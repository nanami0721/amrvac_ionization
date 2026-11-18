!> Fixed-grid overset coupling; local spherical metrics are unchanged.
module mod_yinyang
  use mod_yinyang_staggered, only: yy_staggered_table, yy_staggered_build
  use mod_yinyang_ct_coupling, only: yy_ct_edge_plan, yy_ct_edge_build, yy_ct_edge_exchange
  use mod_global_parameters
  use mod_physics, only: phys_to_primitive, phys_to_conserved, physics_type
  use mod_comm_lib, only: mpistop
  use mod_forest, only: particle_mesh_epoch
  implicit none
  private
  public :: yy_validate, yy_exchange, yy_basis, yy_position, yy_checkpoint, yy_checkpoint_begin, yy_report
  public :: yy_build_staggered_table, yy_prepare_ct_edges, yy_exchange_ct_edges, yy_ct_stage_plan
  type yy_ct_stage_plan
    type(yy_ct_edge_plan) :: edge
    integer,allocatable :: receivers(:,:)
    double precision,allocatable :: packed(:), received(:)
    integer(kind=8) :: mesh_epoch=-1_8
  end type
  integer(kind=8) :: epoch=-1_8
  integer :: nr=0, ns=0, nt=0, raw_nr=0
  integer(kind=8) :: exchange_calls=0_8
  integer, allocatable :: rec(:,:), slot(:,:), request(:,:), sc(:),rc(:),sd(:),rd(:)
  integer, allocatable :: handles(:), statuses(:,:)
  integer, allocatable :: shead(:),snext(:),rhead(:),rnext(:)
  double precision, allocatable :: weights(:,:), rbasis(:,:,:), sbasis(:,:,:)
  double precision, allocatable :: sendbuf(:,:), recvbuf(:,:), scratch(:^D&,:)
  double precision :: pack_seconds=0.d0, wait_seconds=0.d0, unpack_seconds=0.d0, build_seconds=0.d0
contains
  ! Build an explicitly owned table from the real fixed-grid AMRVAC layout.
  ! Caller owns its lifetime and must rebuild after mesh/ownership changes.
  ! This does not activate CT or write state/ghost faces.
  subroutine yy_build_staggered_table(table,s,family,points)
    type(yy_staggered_table),intent(inout) :: table
    type(state),intent(in) :: s(max_blocks)
    integer,intent(in) :: family
    double precision,intent(in) :: points(:,:)
    integer :: nn(3),bn(3)
    integer,allocatable :: directory(:,:)
    double precision,allocatable :: radial(:)
    double precision :: angular(2,2)
    call yy_staggered_mesh(s,bn,nn,radial,angular,directory)
    call yy_staggered_build(table,MPI_COMM_WORLD,yy_patch,family,bn,nn,radial,angular,directory,points)
  end subroutine

  subroutine yy_staggered_mesh(s,bn,nn,rmax,angular,directory)
    type(state),intent(in) :: s(max_blocks)
    integer,intent(out) :: bn(3),nn(3)
    double precision,allocatable,intent(out) :: rmax(:)
    double precision,intent(out) :: angular(2,2)
    integer,allocatable,intent(out) :: directory(:,:)
    {^IFTHREED
    integer :: n,g,p,i,gi
    integer,allocatable :: local(:,:),counts(:),displs(:)
    double precision,allocatable :: lower(:),upper(:),rmin(:)
    double precision :: radius
    if(.not.yinyang_enabled.or.refine_max_level/=1.or.lb_automatic) &
      call mpistop('staggered table requires fixed level-one Yin-Yang layout')
    if(stretched_dim(2).or.stretched_dim(3)) call mpistop('staggered table needs uniform angular mesh')
    nn=[domain_nx1,domain_nx2,domain_nx3];bn=[block_nx1,block_nx2,block_nx3]
    allocate(local(6,igridstail),counts(yy_world_size),displs(yy_world_size))
    allocate(lower(0:nn(1)),upper(0:nn(1)),rmin(0:nn(1)),rmax(0:nn(1)))
    lower=huge(1.d0);upper=-huge(1.d0)
    do n=1,igridstail
      g=igrids(n)
      local(:,n)=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1,g,yy_world_rank,yy_patch]
      do i=ixMlo1-1,ixMhi1
        gi=(node(pig1_,g)-1)*block_nx1+i-ixMlo1+1
        radius=s(g)%x(i,ixMlo2,ixMlo3,1)+0.5d0*s(g)%dx(i,ixMlo2,ixMlo3,1)
        lower(gi)=min(lower(gi),radius);upper(gi)=max(upper(gi),radius)
      end do
    end do
    call MPI_ALLREDUCE(lower,rmin,nn(1)+1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(upper,rmax,nn(1)+1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(any(abs(rmax-rmin)>1.d-12*max(1.d0,abs(rmax)))) &
      call mpistop('staggered table radial meshes disagree or are incomplete')
    call MPI_ALLGATHER(6*igridstail,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    displs(1)=0
    do p=2,yy_world_size
      displs(p)=displs(p-1)+counts(p-1)
    end do
    allocate(directory(6,sum(counts)/6))
    call MPI_ALLGATHERV(local,6*igridstail,MPI_INTEGER,directory,counts,displs,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    angular(:,1)=[xprobmin2,xprobmin3];angular(:,2)=[xprobmax2,xprobmax3]
    }
    {^IFONED call mpistop('staggered Yin-Yang table requires 3D') }
    {^IFTWOD call mpistop('staggered Yin-Yang table requires 3D') }
  end subroutine

  ! Receiver rows are local block ID, native i,j,k and direction. Select only
  ! desired interface edges; uncomputed UCT workspace may never be a donor.
  subroutine yy_prepare_ct_edges(plan,s,receivers)
    type(yy_ct_stage_plan),intent(out) :: plan
    type(state),intent(in) :: s(max_blocks)
    integer,intent(in) :: receivers(:,:)
    {^IFTHREED
    integer :: nn(3),bn(3),n,g,c,idx(3),lo(3),hi(3)
    integer,allocatable :: directory(:,:)
    double precision,allocatable :: radial(:),centers(:,:),widths(:)
    double precision :: angular(2,2),offset(3)
    if(size(receivers,1)/=5) call mpistop('CT edge receiver shape must be (5,n)')
    call yy_staggered_mesh(s,bn,nn,radial,angular,directory)
    allocate(centers(3,size(receivers,2)),widths(size(receivers,2)))
    do n=1,size(receivers,2)
      g=receivers(1,n);c=receivers(5,n);idx=receivers(2:4,n)
      if(c<1.or.c>3) call mpistop('invalid receiving edge direction')
      if(.not.any(igrids(1:igridstail)==g)) call mpistop('receiving CT block is not local')
      lo=[ixMlo1-1,ixMlo2-1,ixMlo3-1];hi=[ixMhi1,ixMhi2,ixMhi3];lo(c)=lo(c)+1
      if(any(idx<lo).or.any(idx>hi)) call mpistop('receiver outside valid UCT stencil')
      offset=0.5d0;offset(c)=0.d0
      centers(:,n)=s(g)%x(idx(1),idx(2),idx(3),:)+offset*s(g)%dx(idx(1),idx(2),idx(3),:)
      widths(n)=s(g)%dx(idx(1),idx(2),idx(3),c)
    end do
    call yy_ct_edge_build(plan%edge,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
      centers,widths,receivers(5,:))
    if(any(.not.plan%edge%covered)) call mpistop('CT edge lacks full donors: boundary closure required')
    do n=1,size(plan%edge%samples%request,2)
      idx=plan%edge%samples%request(2:4,n)+[ixMlo1,ixMlo2,ixMlo3]
      c=plan%edge%samples%request(5,n)
      lo=[ixMlo1-1,ixMlo2-1,ixMlo3-1];lo(c)=lo(c)+1
      if(any(idx<lo).or.any(idx>[ixMhi1,ixMhi2,ixMhi3])) &
        call mpistop('donor outside computed UCT edge stencil')
    end do
    plan%receivers=receivers
    allocate(plan%packed(size(plan%edge%samples%request,2)),plan%received(size(receivers,2)))
    plan%mesh_epoch=particle_mesh_epoch
    }
    {^IFONED call mpistop('CT edge exchange requires 3D') }
    {^IFTWOD call mpistop('CT edge exchange requires 3D') }
  end subroutine

  ! Invoke at the all-block RK-stage boundary, before any edge is overwritten.
  ! All requested donor integrals are packed before the collective returns.
  subroutine yy_exchange_ct_edges(plan,fluxes)
    type(yy_ct_stage_plan),intent(inout) :: plan
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    {^IFTHREED
    integer :: n,g,c,i,j,k
    if(.not.allocated(plan%packed)) call mpistop('CT edge plan has not been built')
    if(size(fluxes)<max_blocks) call mpistop('CT stage block array is too small')
    if(plan%mesh_epoch/=particle_mesh_epoch) call mpistop('CT edge plan is stale')
    do n=1,size(plan%packed)
      g=plan%edge%samples%request(1,n);c=plan%edge%samples%request(5,n)
      i=plan%edge%samples%request(2,n)+ixMlo1
      j=plan%edge%samples%request(3,n)+ixMlo2
      k=plan%edge%samples%request(4,n)+ixMlo3
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT donor stage workspace missing')
      plan%packed(n)=fluxes(g)%edge(i,j,k,c)
    end do
    call yy_ct_edge_exchange(plan%edge,plan%packed,plan%received)
    do n=1,size(plan%received)
      if(.not.plan%edge%covered(n)) cycle
      g=plan%receivers(1,n);i=plan%receivers(2,n);j=plan%receivers(3,n)
      k=plan%receivers(4,n);c=plan%receivers(5,n)
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT receiving workspace missing')
      fluxes(g)%edge(i,j,k,c)=plan%received(n)
    end do
    }
    {^IFONED call mpistop('CT edge exchange requires 3D') }
    {^IFTWOD call mpistop('CT edge exchange requires 3D') }
  end subroutine

  ! Invalidate a previous completion record BEFORE overwriting either snapshot.
  ! Called collectively on the world communicator, outside block/OpenMP loops.
  subroutine yy_checkpoint_begin(filename)
    character(len=*),intent(in) :: filename
    integer :: u,ios
    logical :: exists
    if(.not.yinyang_enabled) return
    if(mype==0) then
      inquire(file=trim(filename)//'.yy.complete',exist=exists)
      if(exists) then
        open(newunit=u,file=trim(filename)//'.yy.complete',status='old',iostat=ios)
        if(ios/=0) call mpistop('cannot invalidate yinyang checkpoint completion')
        close(u,status='delete',iostat=ios)
        if(ios/=0) call mpistop('cannot delete yinyang checkpoint completion')
      end if
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
  end subroutine

  ! Version 2 completion records are published only after both data files and
  ! sidecars are closed. Size checks detect truncation, not same-size corruption.
  subroutine yy_checkpoint(filename,writing)
    character(len=*),intent(in) :: filename
    logical,intent(in) :: writing
    integer :: u,ios,header(8),expected(8),imin,imax,version
    integer(kind=8) :: bytes,saved_bytes
    double precision :: stamp,tmin,tmax
    if(.not.yinyang_enabled) return
    {^IFTHREED
    expected=[2,yy_patch,yy_world_size,domain_nx1,domain_nx2,domain_nx3, &
      merge(1,0,stagger_grid),merge(1,0,B0field)]
    call MPI_ALLREDUCE(it,imin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(it,imax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(global_time,tmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(global_time,tmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(imin/=imax.or.tmin/=tmax) call mpistop('yinyang checkpoint patch times disagree')
    if(mype==0) then
      inquire(file=trim(filename),size=bytes,iostat=ios)
      if(ios/=0.or.bytes<=0) call mpistop('invalid yinyang checkpoint data size')
      if(writing) then
        open(newunit=u,file=trim(filename)//'.yy',status='replace',iostat=ios)
        if(ios/=0) call mpistop('cannot write yinyang checkpoint metadata')
        write(u,*,iostat=ios) expected
        if(ios/=0) call mpistop('cannot write yinyang checkpoint header')
        write(u,*,iostat=ios) it,global_time,bytes
        if(ios/=0) call mpistop('cannot write yinyang checkpoint time and size')
      else
        open(newunit=u,file=trim(filename)//'.yy',status='old',iostat=ios)
        if(ios/=0) call mpistop('missing yinyang checkpoint metadata')
        read(u,*,iostat=ios) version
        if(ios/=0) call mpistop('invalid yinyang checkpoint metadata')
        rewind(u)
        if(version==1) then
          read(u,*,iostat=ios) header(1:6)
          if(ios/=0) call mpistop('invalid legacy yinyang checkpoint metadata')
          if(any(header(2:6)/=expected(2:6))) call mpistop('legacy yinyang checkpoint mesh or MPI mismatch')
          read(u,*,iostat=ios) imin,stamp
          if(ios/=0) call mpistop('invalid legacy yinyang checkpoint time')
          if(yy_world_rank==0) write(*,'(a)') &
            'Yin-Yang: legacy v1 checkpoint has no completion or file-size guarantee'
        else if(version==2) then
          read(u,*,iostat=ios) header
          if(ios/=0) call mpistop('invalid yinyang checkpoint header')
          if(any(header/=expected)) call mpistop('yinyang checkpoint mesh, patch, MPI or field mode mismatch')
          read(u,*,iostat=ios) imin,stamp,saved_bytes
          if(ios/=0) call mpistop('invalid yinyang checkpoint time and size')
          if(saved_bytes/=bytes) call mpistop('yinyang checkpoint data size mismatch')
        else
          call mpistop('unsupported yinyang checkpoint version')
        end if
        if(imin/=it.or.stamp/=global_time) call mpistop('yinyang checkpoint timestamp mismatch')
      end if
      close(u,iostat=ios)
      if(ios/=0) call mpistop('cannot close yinyang checkpoint metadata')
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
    if(mype==0) then
      if(writing) then
        open(newunit=u,file=trim(filename)//'.yy.complete',status='replace',iostat=ios)
        if(ios/=0) call mpistop('cannot publish yinyang checkpoint completion')
        write(u,*,iostat=ios) it,global_time,bytes
        if(ios/=0) call mpistop('cannot write yinyang checkpoint completion')
        close(u,iostat=ios)
        if(ios/=0) call mpistop('cannot close yinyang checkpoint completion')
      else if(version==2) then
        open(newunit=u,file=trim(filename)//'.yy.complete',status='old',iostat=ios)
        if(ios/=0) call mpistop('missing yinyang checkpoint completion')
        read(u,*,iostat=ios) imin,stamp,saved_bytes
        if(ios/=0) call mpistop('invalid yinyang checkpoint completion')
        close(u)
        if(imin/=it.or.stamp/=global_time.or.saved_bytes/=bytes) &
          call mpistop('stale yinyang checkpoint completion')
      end if
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
    }
  end subroutine

  subroutine yy_report()
    double precision :: local(4),maximum(4),average(4)
    integer(kind=8) :: samples(2),totals(2),calls_min,calls_max
    if(.not.yinyang_enabled) return
    local=[build_seconds,pack_seconds,wait_seconds,unpack_seconds]
    call MPI_REDUCE(local,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(local,average,4,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierrmpi)
    samples=[int(raw_nr,8),int(nr,8)]
    call MPI_REDUCE(samples,totals,2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(exchange_calls,calls_min,1,MPI_INTEGER8,MPI_MIN,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(exchange_calls,calls_max,1,MPI_INTEGER8,MPI_MAX,0,MPI_COMM_WORLD,ierrmpi)
    if(yy_world_rank==0) then
      average=average/yy_world_size
      write(*,'(a,4es14.6)') 'Yin-Yang max-rank build / pack / exchange / unpack seconds: ',maximum
      write(*,'(a,4es14.6)') 'Yin-Yang mean-rank build / pack / exchange / unpack seconds: ',average
      write(*,'(a,2i14)') 'Yin-Yang donor samples raw / unique per exchange: ',totals
      write(*,'(a,2i14)') 'Yin-Yang payload bytes raw / unique per exchange: ',totals*int(nw,8)*8_8
      write(*,'(a,2i14)') 'Yin-Yang exchange calls min / max: ',calls_min,calls_max
    end if
  end subroutine

  ! Deduplicate requests for each sending peer. Only schedule construction uses
  ! hashing; runtime slot lookups and their four-term summation order are unchanged.
  subroutine yy_unique_requests(outreq)
    integer,allocatable,intent(inout) :: outreq(:,:)
    integer,allocatable :: compact(:,:),remap(:),table(:),peer(:),unique_counts(:)
    integer :: capacity,p,j,q,h,found,nunique,t
    integer(kind=8) :: key
    raw_nr=nr
    if(nr>(huge(nr)-1)/2) call mpistop('yinyang request table exceeds integer capacity')
    capacity=max(1,2*nr+1)
    allocate(compact(4,nr),remap(nr),table(capacity),peer(nr),unique_counts(yy_world_size))
    table=0;unique_counts=0;nunique=0
    do p=1,yy_world_size
      do j=sd(p)+1,sd(p)+sc(p)
        key=int(p,8)
        do q=1,4
          key=modulo(key*1000003_8+int(outreq(q,j),8),int(capacity,8))
        end do
        h=int(key)+1
        do
          found=table(h)
          if(found==0) then
            nunique=nunique+1
            compact(:,nunique)=outreq(:,j);peer(nunique)=p
            table(h)=nunique;remap(j)=nunique
            unique_counts(p)=unique_counts(p)+1
            exit
          end if
          if(peer(found)==p) then
            if(all(compact(:,found)==outreq(:,j))) then
              remap(j)=found
              exit
            end if
          end if
          h=mod(h,capacity)+1
        end do
      end do
    end do
    do t=1,nt
      do q=1,4
        slot(q,t)=remap(slot(q,t))
      end do
    end do
    sc=unique_counts;sd(1)=0
    do p=2,yy_world_size
      sd(p)=sd(p-1)+sc(p-1)
    end do
    nr=nunique
    deallocate(outreq)
    allocate(outreq(4,nr))
    outreq=compact(:,1:nr)
  end subroutine

  pure subroutine yy_basis(theta,phi,patch,b)
    double precision,intent(in) :: theta,phi
    integer,intent(in) :: patch
    double precision,intent(out) :: b(3,3)
    double precision :: tmp(3)
    b(:,1)=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
    b(:,2)=[cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)]
    b(:,3)=[-sin(phi),cos(phi),0.d0]
    if(patch==1) then
      b(1,:)=-b(1,:)
      tmp=b(2,:); b(2,:)=b(3,:); b(3,:)=tmp
    end if
  end subroutine

  pure subroutine yy_position(r,theta,phi,patch,v)
    double precision,intent(in) :: r,theta,phi
    integer,intent(in) :: patch
    double precision,intent(out) :: v(3)
    double precision :: b(3,3)
    call yy_basis(theta,phi,patch,b)
    v=r*b(:,1)
  end subroutine

  subroutine yy_validate()
    if(.not.yinyang_enabled) return
    {^IFTHREED
    if(refine_max_level/=1.or.lb_automatic.or.local_timestep) &
      call mpistop('yinyang v1 requires fixed level-one mesh and global time stepping')
    if(stagger_grid.or.B0field.or.use_particles.or.use_multigrid) &
      call mpistop('yinyang v1: CT, B0, particles and multigrid are unsupported')
    if(any(periodB).or.any(aperiodB)) call mpistop('yinyang has no periodic boundary')
    if(xprobmin1<=0.d0.or.xprobmin2<=0.d0.or.xprobmax2>=dpi) &
      call mpistop('yinyang requires a positive-radius shell excluding polar axes')
    if(stretched_dim(2).or.stretched_dim(3)) call mpistop('yinyang angular mesh must be uniform')
    if(physics_type/='mhd'.and.physics_type/='hd'.and.physics_type/='rho') &
      call mpistop('yinyang supports rho, hd and ideal GLM mhd only')
    if(nw/=nwgc.or.iwstart/=1.or.nwaux>0.or.nwextra>0) &
      call mpistop('yinyang requires complete primary-state exchange')
    if(yy_world_rank==0) write(*,'(a)') 'Yin-Yang v1: interpolation coupling (not flux conservative)'
    }
  end subroutine

  subroutine yy_build(s)
    type(state),intent(in) :: s(max_blocks)
    {^IFTHREED
    integer :: nblocks,nlocal,n,g,a,b,j,k,i,q,t,p,off,dj,dk,bi,bj,bk,ci,cj,ck,pass,gi,gj,gk
    integer,allocatable :: counts(:),displs(:),local(:,:),allblocks(:,:),owners(:,:,:,:),ids(:,:,:,:)
    integer,allocatable :: outreq(:,:),cursor(:),sc4(:),rc4(:),sd4(:),rd4(:)
    double precision :: v(3),theta,phi,aa,bb,wa,wb,dteta,dphi,mat(3,3)
    if(allocated(rec)) then
      deallocate(rec,slot,request,sc,rc,sd,rd,handles,statuses,weights,rbasis,sbasis,sendbuf,recvbuf,scratch,shead,snext,rhead,rnext)
    end if
    allocate(counts(yy_world_size),displs(yy_world_size))
    nlocal=igridstail
    allocate(local(5,nlocal))
    do n=1,nlocal
      g=igrids(n)
      local(:,n)=[node(pig1_,g),node(pig2_,g),node(pig3_,g),g,yy_world_rank]
    end do
    call MPI_ALLGATHER(5*nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    displs(1)=0
    do p=2,yy_world_size
      displs(p)=displs(p-1)+counts(p-1)
    end do
    nblocks=sum(counts)/5
    allocate(allblocks(5,nblocks))
    call MPI_ALLGATHERV(local,5*nlocal,MPI_INTEGER,allblocks,counts,displs,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    allocate(owners(ng1(1),ng2(1),ng3(1),0:1),ids(ng1(1),ng2(1),ng3(1),0:1))
    owners=-1; ids=-1
    do n=1,nblocks
      p=allblocks(5,n)/(yy_world_size/2)
      bi=allblocks(1,n);bj=allblocks(2,n);bk=allblocks(3,n)
      owners(bi,bj,bk,p)=allblocks(5,n)
      ids(bi,bj,bk,p)=allblocks(4,n)
    end do
    if(any(owners<0)) call mpistop('yinyang missing level-one block')
    allocate(sc(yy_world_size),rc(yy_world_size),sd(yy_world_size),rd(yy_world_size),cursor(yy_world_size))
    sc=0; nt=0
    dteta=(xprobmax2-xprobmin2)/domain_nx2
    dphi=(xprobmax3-xprobmin3)/domain_nx3
    do pass=1,2
      t=0
      do n=1,nlocal
        g=igrids(n)
        do k=ixGlo3,ixGhi3
          gk=(node(pig3_,g)-1)*block_nx3+k-ixMlo3
          do j=ixGlo2,ixGhi2
            gj=(node(pig2_,g)-1)*block_nx2+j-ixMlo2
            if(gj>=0.and.gj<domain_nx2.and.gk>=0.and.gk<domain_nx3) cycle
            theta=s(g)%x(ixMlo1,j,k,2); phi=s(g)%x(ixMlo1,j,k,3)
            v=[-sin(theta)*cos(phi),cos(theta),sin(theta)*sin(phi)]
            theta=atan2(sqrt(v(1)**2+v(2)**2),v(3));phi=atan2(v(2),v(1))
            aa=(theta-xprobmin2)/dteta-0.5d0;bb=(phi-xprobmin3)/dphi-0.5d0
            dj=floor(aa);dk=floor(bb);wa=aa-dj;wb=bb-dk
            if(dj<0.or.dj+1>=domain_nx2.or.dk<0.or.dk+1>=domain_nx3) &
              call mpistop('yinyang ghost lacks interior donors: increase angular overlap')
            do i=ixGlo1,ixGhi1
              gi=(node(pig1_,g)-1)*block_nx1+i-ixMlo1
              bi=max(1,min(ng1(1),gi/block_nx1+1));ci=gi-(bi-1)*block_nx1+ixMlo1
              t=t+1
              if(pass==2) then
                rec(:,t)=[g,i,j,k]
                call yy_basis(s(g)%x(i,j,k,2),s(g)%x(i,j,k,3),yy_patch,mat)
                rbasis(:,:,t)=transpose(mat)
              end if
              q=0
              do b=0,1
                do a=0,1
                  q=q+1
                  bj=(dj+a)/block_nx2+1;bk=(dk+b)/block_nx3+1
                  cj=mod(dj+a,block_nx2)+ixMlo2;ck=mod(dk+b,block_nx3)+ixMlo3
                  p=owners(bi,bj,bk,1-yy_patch)+1
                  if(pass==1) then
                    sc(p)=sc(p)+1
                  else
                    off=sd(p)+cursor(p)+1;cursor(p)=cursor(p)+1
                    slot(q,t)=off
                    weights(q,t)=merge(1.d0-wa,wa,a==0)*merge(1.d0-wb,wb,b==0)
                    outreq(:,off)=[ids(bi,bj,bk,1-yy_patch),ci,cj,ck]
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
      if(pass==1) then
        nt=t;nr=sum(sc)
        allocate(rec(4,nt),slot(4,nt),weights(4,nt),rbasis(3,3,nt),outreq(4,nr))
        sd(1)=0
        do p=2,yy_world_size
          sd(p)=sd(p-1)+sc(p-1)
        end do
        cursor=0
      end if
    end do
    call yy_unique_requests(outreq)
    call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    rd(1)=0
    do p=2,yy_world_size
      rd(p)=rd(p-1)+rc(p-1)
    end do
    ns=sum(rc)
    allocate(request(4,ns),sbasis(3,3,ns),sendbuf(nw,ns),recvbuf(nw,nr))
    allocate(sc4(yy_world_size),rc4(yy_world_size),sd4(yy_world_size),rd4(yy_world_size))
    sc4=4*sc;rc4=4*rc;sd4=4*sd;rd4=4*rd
    call MPI_ALLTOALLV(outreq,sc4,sd4,MPI_INTEGER,request,rc4,rd4,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    do t=1,ns
      g=request(1,t);i=request(2,t);j=request(3,t);k=request(4,t)
      call yy_basis(s(g)%x(i,j,k,2),s(g)%x(i,j,k,3),yy_patch,sbasis(:,:,t))
    end do
    allocate(scratch(ixG^T,nw),handles(2*yy_world_size),statuses(MPI_STATUS_SIZE,2*yy_world_size))
    allocate(shead(max_blocks),rhead(max_blocks),snext(ns),rnext(nt))
    shead=0;rhead=0
    do t=ns,1,-1
      g=request(1,t);snext(t)=shead(g);shead(g)=t
    end do
    do t=nt,1,-1
      g=rec(1,t);rnext(t)=rhead(g);rhead(g)=t
    end do
    epoch=particle_mesh_epoch
    if(yy_world_rank==0) write(*,'(a,i0,a,i0)') 'Yin-Yang schedule: rank 0 receivers=',nt,', samples=',ns
    }
  end subroutine

  subroutine yy_exchange(s,nwstart,nwbc)
    type(state),target,intent(inout) :: s(max_blocks)
    integer,intent(in) :: nwstart,nwbc
    {^IFTHREED
    integer :: n,g,t,i,j,k,q,iv,iw,p,nh,side,axis,ixO^L,lo(3),hi(3),bn(3),dn(3)
    double precision :: v(3),val(nw),t0
    if(nwstart/=1.or.nwbc/=nw) call mpistop('yinyang partial-state exchange unsupported')
    if(epoch/=particle_mesh_epoch) then
      t0=MPI_WTIME()
      call yy_build(s)
      build_seconds=build_seconds+MPI_WTIME()-t0
    end if
    exchange_calls=exchange_calls+1_8
    bn=[block_nx1,block_nx2,block_nx3];dn=[domain_nx1,domain_nx2,domain_nx3]
    t0=MPI_WTIME()
    do n=1,igridstail
      g=igrids(n)
      if(shead(g)==0) cycle
      block=>s(g)
      scratch=s(g)%w
      ixO^L=ixG^LL;
      ixOmin2=max(ixGlo2,ixMlo2-(node(pig2_,g)-1)*block_nx2)
      ixOmax2=min(ixGhi2,ixMlo2+domain_nx2-1-(node(pig2_,g)-1)*block_nx2)
      ixOmin3=max(ixGlo3,ixMlo3-(node(pig3_,g)-1)*block_nx3)
      ixOmax3=min(ixGhi3,ixMlo3+domain_nx3-1-(node(pig3_,g)-1)*block_nx3)
      call phys_to_primitive(ixG^LL,ixO^L,scratch,s(g)%x)
      t=shead(g)
      do while(t/=0)
        i=request(2,t);j=request(3,t);k=request(4,t)
        sendbuf(:,t)=scratch(i,j,k,:)
        do iv=1,nvector
          iw=iw_vector(iv);v=scratch(i,j,k,iw+1:iw+3)
          sendbuf(iw+1:iw+3,t)=matmul(sbasis(:,:,t),v)
        end do
        t=snext(t)
      end do
    end do
    pack_seconds=pack_seconds+MPI_WTIME()-t0
    nh=0;t0=MPI_WTIME()
    do p=1,yy_world_size
      if(sc(p)>0) then
        nh=nh+1
        call MPI_IRECV(recvbuf(1,sd(p)+1),nw*sc(p),MPI_DOUBLE_PRECISION,p-1,27131,MPI_COMM_WORLD,handles(nh),ierrmpi)
      end if
      if(rc(p)>0) then
        nh=nh+1
        call MPI_ISEND(sendbuf(1,rd(p)+1),nw*rc(p),MPI_DOUBLE_PRECISION,p-1,27131,MPI_COMM_WORLD,handles(nh),ierrmpi)
      end if
    end do
    if(nh>0) call MPI_WAITALL(nh,handles,statuses,ierrmpi)
    wait_seconds=wait_seconds+MPI_WTIME()-t0
    t0=MPI_WTIME()
    do n=1,igridstail
      g=igrids(n)
      if(rhead(g)==0) cycle
      block=>s(g)
      t=rhead(g)
      do while(t/=0)
        i=rec(2,t);j=rec(3,t);k=rec(4,t)
        val=0.d0
        do q=1,4
          val=val+weights(q,t)*recvbuf(:,slot(q,t))
        end do
        do iv=1,nvector
          iw=iw_vector(iv);v=val(iw+1:iw+3)
          val(iw+1:iw+3)=matmul(rbasis(:,:,t),v)
        end do
        scratch(i,j,k,:)=val
        t=rnext(t)
      end do
      ! Convert rectangular ghost slabs in bulk; every read cell was filled above.
      do axis=2,3
        do side=1,2
          lo=[ixGlo1,ixGlo2,ixGlo3];hi=[ixGhi1,ixGhi2,ixGhi3]
          if(side==1) then
            hi(axis)=nghostcells-(node(plevel_+axis,g)-1)*bn(axis)
          else
            lo(axis)=nghostcells+dn(axis)-(node(plevel_+axis,g)-1)*bn(axis)+1
          end if
          if(axis==3) then
            lo(2)=max(ixGlo2,ixMlo2-(node(pig2_,g)-1)*block_nx2)
            hi(2)=min(ixGhi2,ixMlo2+domain_nx2-1-(node(pig2_,g)-1)*block_nx2)
          end if
          if(any(lo>hi)) cycle
          ixOmin1=lo(1);ixOmin2=lo(2);ixOmin3=lo(3)
          ixOmax1=hi(1);ixOmax2=hi(2);ixOmax3=hi(3)
          call phys_to_conserved(ixG^LL,ixO^L,scratch,s(g)%x)
          s(g)%w(ixO^S,:)=scratch(ixO^S,:)
        end do
      end do
    end do
    unpack_seconds=unpack_seconds+MPI_WTIME()-t0
    }
  end subroutine
end module mod_yinyang
