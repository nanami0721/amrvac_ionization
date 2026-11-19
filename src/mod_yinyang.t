!> Fixed-grid overset coupling; local spherical metrics are unchanged.
module mod_yinyang
  use mod_yinyang_staggered, only: yy_staggered_table, yy_staggered_build
  use mod_yinyang_ct_coupling, only: yy_ct_edge_plan, yy_ct_edge_build, yy_ct_edge_exchange, &
    yy_ct_face_plan, yy_ct_face_build, yy_ct_face_exchange, yy_ct_close_ghosts
  use mod_yinyang_patch_sync, only: yy_patch_copy_plan, yy_patch_copy_build, yy_patch_copy_exchange, &
    yy_sync_edges, yy_sync_faces, yy_sync_extended_faces
  use mod_global_parameters
  use mod_physics, only: phys_to_primitive, phys_to_conserved, physics_type, &
    phys_ct_sync_center,phys_get_ct_edges,phys_apply_ct_edges
  use mod_comm_lib, only: mpistop
  use mod_forest, only: particle_mesh_epoch
  implicit none
  private
  ! Owned here: only input parsing and Yin-Yang initialization use this option.
  logical,public :: yy_ct_face_moments=.false.
  public :: yy_validate, yy_exchange, yy_basis, yy_position, yy_checkpoint, yy_checkpoint_begin, yy_report
  public :: yy_build_staggered_table, yy_prepare_ct_edges, yy_exchange_ct_edges, yy_ct_stage_plan
  public :: yy_ct_ghost_plan, yy_prepare_ct_ghosts, yy_fill_ct_ghosts
  type yy_ct_ghost_plan
    type(yy_ct_face_plan) :: faces
    integer,allocatable :: receivers(:,:)
    double precision,allocatable :: radial(:),tangent(:),values(:)
    integer(kind=8) :: mesh_epoch=-1_8
  end type
  public :: yy_ct_patch_plan, yy_prepare_ct_patch, yy_sync_ct_edges, yy_ordered_ct_edges
  public :: yy_sync_ct_faces, yy_ordered_ct_ghosts, yy_ct_edge_is_owner
  public :: yy_complete_ct_faces, yy_sync_ct_state, yy_ct_radial_faces
  type yy_ct_patch_plan
    type(yy_patch_copy_plan) :: edge,interior_faces,extended_faces
    integer(kind=8) :: mesh_epoch=-1_8
  end type
  type yy_ct_stage_plan
    type(yy_ct_edge_plan) :: edge
    integer,allocatable :: receivers(:,:)
    double precision,allocatable :: packed(:), received(:)
    integer(kind=8) :: mesh_epoch=-1_8
  end type
  public :: yy_runtime_ct_edges, yy_runtime_ct_ghosts, yy_ct_radial_cells
  type(yy_ct_patch_plan),save :: runtime_patch
  type(yy_ct_stage_plan),save :: runtime_edges
  type(yy_ct_ghost_plan),save :: runtime_ghosts
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
  subroutine yy_prepare_ct_edges(plan,s,receivers,radial_extrapolate,native_midpoint)
    type(yy_ct_stage_plan),intent(out) :: plan
    type(state),intent(in) :: s(max_blocks)
    integer,intent(in) :: receivers(:,:)
    logical,optional,intent(in) :: radial_extrapolate,native_midpoint
    {^IFTHREED
    logical :: midpoint
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
    midpoint=.true.;if(present(native_midpoint))midpoint=native_midpoint
    call yy_ct_edge_build(plan%edge,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
      centers,widths,receivers(5,:),radial_extrapolate,native_midpoint=midpoint)
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

  ! Enumerate the union of faces bounding external angular ghost cells. Active
  ! faces, including active-boundary anchors, are excluded before interpolation.
  subroutine yy_prepare_ct_ghosts(plan,s,face_moments)
    type(yy_ct_ghost_plan),intent(out) :: plan
    type(state),intent(in) :: s(max_blocks)
    logical,optional,intent(in) :: face_moments
    {^IFTHREED
    integer :: bn(3),nn(3),n,g,i,j,k,c,d,side,pass,t,idx(3),above(3),lo(3),hi(3),fl(3),pig(3),nb(3)
    integer,allocatable :: directory(:,:)
    logical :: moments
    logical :: sides(2,2),selected,mask(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3)
    double precision,allocatable :: radial(:),centers(:,:),widths(:,:)
    double precision :: angular(2,2),offset(3)
    call yy_staggered_mesh(s,bn,nn,radial,angular,directory)
    lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3];nb=nn/bn
    do pass=1,2
      t=0
      do n=1,igridstail
        g=igrids(n);pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
        sides(1,:)=(pig(2:3)==1);sides(2,:)=(pig(2:3)==nb(2:3));mask=.false.
        do k=ixGlo3,ixGhi3
          do j=ixGlo2,ixGhi2
            selected=((j<lo(2).and.sides(1,1)).or.(j>hi(2).and.sides(2,1))) &
              .and.k>=lo(3).and.k<=hi(3)
            if((k<lo(3).and.sides(1,2)).or.(k>hi(3).and.sides(2,2))) &
              selected=(j>=lo(2).and.j<=hi(2)).or.(j<lo(2).and.sides(1,1)).or.(j>hi(2).and.sides(2,1))
            if(.not.selected) cycle
            do i=lo(1),hi(1)
              do c=1,3
                do side=0,1
                  idx=[i,j,k];idx(c)=idx(c)-side
                  fl=lo;fl(c)=fl(c)-1
                  if(all(idx>=fl).and.all(idx<=hi)) cycle
                  mask(idx(1),idx(2),idx(3),c)=.true.
                end do
              end do
            end do
          end do
        end do
        do c=1,3
          do k=ixGlo3-1,ixGhi3
            do j=ixGlo2-1,ixGhi2
              do i=ixGlo1-1,ixGhi1
                if(.not.mask(i,j,k,c)) cycle
                t=t+1
                if(pass==1) cycle
                plan%receivers(:,t)=[g,i,j,k,c];idx=[i,j,k]
                above=max(idx,[ixGlo1,ixGlo2,ixGlo3])
                offset=0.d0;offset(c)=0.5d0
                ! Lower outer face has no center x at index zero. Locate it
                ! from the adjacent cell center, not from uninitialized x.
                if(idx(c)<above(c)) offset(c)=-0.5d0
                widths(:,t)=s(g)%dx(above(1),above(2),above(3),:)
                centers(:,t)=s(g)%x(above(1),above(2),above(3),:)+ &
                  offset*s(g)%dx(above(1),above(2),above(3),:)
              end do
            end do
          end do
        end do
      end do
      if(pass==1) allocate(plan%receivers(5,t),centers(3,t),widths(3,t),plan%values(t))
    end do
    moments=.false.;if(present(face_moments))moments=face_moments
    if(moments)then
      call yy_ct_face_build(plan%faces,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
        centers,plan%receivers(5,:),widths)
    else
      call yy_ct_face_build(plan%faces,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
        centers,plan%receivers(5,:))
    end if
    if(any(.not.plan%faces%covered)) call mpistop('magnetic ghost face lacks donors: increase overlap')
    allocate(plan%radial(size(plan%faces%radial%request,2)),plan%tangent(size(plan%faces%tangent%request,2)))
    plan%mesh_epoch=particle_mesh_epoch
    }
    {^IFONED call mpistop('magnetic ghost fill requires 3D') }
    {^IFTWOD call mpistop('magnetic ghost fill requires 3D') }
  end subroutine

  ! Fill face B in angular halos, optionally followed by outward flux closure.
  ! Both donor families are packed before any writes, so a ghost never becomes
  ! a donor through ordering. Caller handles physical radial and same-patch BCs.
  subroutine yy_fill_ct_ghosts(plan,s,constrain)
    type(yy_ct_ghost_plan),intent(inout) :: plan
    type(state),intent(inout) :: s(max_blocks)
    logical,intent(in) :: constrain
    {^IFTHREED
    integer :: n,g,i,j,k,c,idx(3),lo(3),hi(3),status,pig(3),nb(3)
    logical :: sides(2,2)
    if(.not.allocated(plan%values)) call mpistop('magnetic ghost plan has not been built')
    if(plan%mesh_epoch/=particle_mesh_epoch) call mpistop('magnetic ghost plan is stale')
    lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3];nb=[ng1(1),ng2(1),ng3(1)]
    do n=1,igridstail
      g=igrids(n)
      if(.not.allocated(s(g)%ws)) call mpistop('magnetic face workspace missing')
    end do
    do n=1,size(plan%radial)
      g=plan%faces%radial%request(1,n);idx=plan%faces%radial%request(2:4,n)+lo
      c=plan%faces%radial%request(5,n)
      plan%radial(n)=s(g)%ws(idx(1),idx(2),idx(3),c)
    end do
    do n=1,size(plan%tangent)
      g=plan%faces%tangent%request(1,n);idx=plan%faces%tangent%request(2:4,n)+lo
      c=plan%faces%tangent%request(5,n)
      plan%tangent(n)=s(g)%ws(idx(1),idx(2),idx(3),c)
    end do
    call yy_ct_face_exchange(plan%faces,plan%radial,plan%tangent,plan%values)
    do n=1,size(plan%values)
      g=plan%receivers(1,n);i=plan%receivers(2,n);j=plan%receivers(3,n)
      k=plan%receivers(4,n);c=plan%receivers(5,n)
      s(g)%ws(i,j,k,c)=plan%values(n)
    end do
    if(constrain) then
      do n=1,igridstail
        g=igrids(n);pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
        sides(1,:)=(pig(2:3)==1);sides(2,:)=(pig(2:3)==nb(2:3))
        call yy_ct_close_ghosts([ixGlo1-1,ixGlo2-1,ixGlo3-1],lo,hi, &
          [ixGlo1,ixGlo2,ixGlo3],[ixGhi1,ixGhi2,ixGhi3],sides,s(g)%ws,s(g)%surfaceC,status)
        if(status/=0) call mpistop('invalid geometry or field in magnetic ghost closure')
      end do
    end if
    }
    {^IFONED call mpistop('magnetic ghost fill requires 3D') }
    {^IFTWOD call mpistop('magnetic ghost fill requires 3D') }
  end subroutine

  subroutine yy_prepare_ct_patch(plan,s)
    type(yy_ct_patch_plan),intent(out) :: plan
    type(state),intent(in) :: s(max_blocks)
    integer :: bn(3),nn(3)
    integer,allocatable :: directory(:,:)
    double precision,allocatable :: radial(:)
    double precision :: angular(2,2)
    call yy_staggered_mesh(s,bn,nn,radial,angular,directory)
    call yy_patch_copy_build(plan%edge,MPI_COMM_WORLD,yy_patch,yy_sync_edges,bn,nn,nghostcells,directory)
    call yy_patch_copy_build(plan%interior_faces,MPI_COMM_WORLD,yy_patch,yy_sync_faces,bn,nn,nghostcells,directory)
    call yy_patch_copy_build(plan%extended_faces,MPI_COMM_WORLD,yy_patch,yy_sync_extended_faces, &
      bn,nn,nghostcells,directory)
    plan%mesh_epoch=particle_mesh_epoch
  end subroutine

  logical function yy_ct_edge_is_owner(g,idx,c) result(owner)
    integer,intent(in) :: g,idx(3),c
    {^IFTHREED
    integer :: d,lo(3),pig(3)
    lo=[ixMlo1,ixMlo2,ixMlo3];pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
    owner=.true.
    do d=1,3
      if(d==c) cycle
      if(idx(d)==lo(d)-1.and.pig(d)>1) owner=.false.
    end do
    }
    {^IFONED owner=.false. }
    {^IFTWOD owner=.false. }
  end function

  subroutine yy_sync_ct_edges(plan,fluxes)
    type(yy_ct_patch_plan),intent(inout) :: plan
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    {^IFTHREED
    integer :: n,g,c,idx(3),lo(3)
    if(plan%mesh_epoch/=particle_mesh_epoch) call mpistop('CT patch plan is stale')
    if(.not.allocated(plan%edge%send)) call mpistop('CT patch plan has not been built')
    if(size(fluxes)<max_blocks) call mpistop('CT stage array is too small')
    lo=[ixMlo1,ixMlo2,ixMlo3]
    do n=1,size(plan%edge%send)
      g=plan%edge%request(1,n);idx=plan%edge%request(2:4,n)+lo;c=plan%edge%request(5,n)
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT patch donor workspace missing')
      plan%edge%send(n)=fluxes(g)%edge(idx(1),idx(2),idx(3),c)
    end do
    call yy_patch_copy_exchange(plan%edge)
    do n=1,size(plan%edge%recv)
      g=plan%edge%receivers(1,n);idx=plan%edge%receivers(2:4,n)+lo;c=plan%edge%receivers(5,n)
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT patch receiver workspace missing')
      fluxes(g)%edge(idx(1),idx(2),idx(3),c)=plan%edge%recv(n)
    end do
    }
    {^IFONED call mpistop('CT patch sync requires 3D') }
    {^IFTWOD call mpistop('CT patch sync requires 3D') }
  end subroutine

  ! The overset receiver list must select geometric owners. Otherwise the final
  ! same-patch broadcast would erase an update applied only to a non-owner copy.
  subroutine yy_ordered_ct_edges(patch_plan,cross_plan,fluxes)
    type(yy_ct_patch_plan),intent(inout) :: patch_plan
    type(yy_ct_stage_plan),intent(inout) :: cross_plan
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    integer :: n
    if(.not.allocated(cross_plan%receivers)) call mpistop('CT cross-patch plan has not been built')
    do n=1,size(cross_plan%receivers,2)
      if(.not.yy_ct_edge_is_owner(cross_plan%receivers(1,n),cross_plan%receivers(2:4,n), &
        cross_plan%receivers(5,n))) call mpistop('CT overset receiver must be the geometric owner')
    end do
    call yy_sync_ct_edges(patch_plan,fluxes)
    call yy_exchange_ct_edges(cross_plan,fluxes)
    call yy_sync_ct_edges(patch_plan,fluxes)
  end subroutine

  subroutine yy_sync_ct_faces(plan,s,extended)
    type(yy_ct_patch_plan),target,intent(inout) :: plan
    type(state),intent(inout) :: s(max_blocks)
    logical,intent(in) :: extended
    {^IFTHREED
    type(yy_patch_copy_plan),pointer :: table
    integer :: n,g,c,idx(3),lo(3)
    if(plan%mesh_epoch/=particle_mesh_epoch) call mpistop('CT patch plan is stale')
    if(extended) then
      table=>plan%extended_faces
    else
      table=>plan%interior_faces
    end if
    if(.not.allocated(table%send)) call mpistop('CT face sync plan has not been built')
    lo=[ixMlo1,ixMlo2,ixMlo3]
    do n=1,size(table%send)
      g=table%request(1,n);idx=table%request(2:4,n)+lo;c=table%request(5,n)
      if(.not.allocated(s(g)%ws)) call mpistop('CT face donor state missing')
      table%send(n)=s(g)%ws(idx(1),idx(2),idx(3),c)
    end do
    call yy_patch_copy_exchange(table)
    do n=1,size(table%recv)
      g=table%receivers(1,n);idx=table%receivers(2:4,n)+lo;c=table%receivers(5,n)
      if(.not.allocated(s(g)%ws)) call mpistop('CT face receiver state missing')
      s(g)%ws(idx(1),idx(2),idx(3),c)=table%recv(n)
    end do
    }
    {^IFONED call mpistop('CT face sync requires 3D') }
    {^IFTWOD call mpistop('CT face sync requires 3D') }
  end subroutine

  subroutine yy_ordered_ct_ghosts(patch_plan,cross_plan,s)
    type(yy_ct_patch_plan),intent(inout) :: patch_plan
    type(yy_ct_ghost_plan),intent(inout) :: cross_plan
    type(state),intent(inout) :: s(max_blocks)
    call yy_sync_ct_faces(patch_plan,s,.false.)
    call yy_fill_ct_ghosts(cross_plan,s,.true.)
    call yy_sync_ct_faces(patch_plan,s,.true.)
  end subroutine

  ! First finish angular faces, then radial physical halos across the entire
  ! angular extent, including corners. Cell thermodynamic halos must already be
  ! populated (ordinary/overset exchange) before this explicit routine is called.
  subroutine yy_complete_ct_faces(patch_plan,cross_plan,s)
    type(yy_ct_patch_plan),intent(inout) :: patch_plan
    type(yy_ct_ghost_plan),intent(inout) :: cross_plan
    type(state),intent(inout) :: s(max_blocks)
    call yy_ordered_ct_ghosts(patch_plan,cross_plan,s)
    call yy_ct_radial_faces(s)
  end subroutine

  ! Explicit radial face policy: constant tangential B,
  ! zero finite-volume divergence in radial ghost cells. Boundary Br is anchored.
  ! This is not a characteristic outflow or a line-tied velocity boundary.
  subroutine yy_ct_radial_faces(s)
    type(state),intent(inout) :: s(max_blocks)
    {^IFTHREED
    integer :: n,g,side,anchor,first,last,stride,i,j,k,c,jlo,klo
    double precision :: transverse
    if(.not.yinyang_enabled.or.refine_max_level/=1.or.lb_automatic) &
      call mpistop('CT radial fill requires fixed Yin-Yang mesh')
    do n=1,igridstail
      g=igrids(n)
      if(.not.allocated(s(g)%ws)) call mpistop('CT radial face state missing')
      do side=1,2
        if(side==1.and.node(pig1_,g)/=1) cycle
        if(side==2.and.node(pig1_,g)/=domain_nx1/block_nx1) cycle
        anchor=ixMlo1;first=ixMlo1-1;last=ixGlo1;stride=-1
        if(side==2) then
          anchor=ixMhi1;first=ixMhi1+1;last=ixGhi1;stride=1
        end if
        ! Angular face halos at the anchor radius are already complete.
        do c=2,3
          jlo=ixGlo2;klo=ixGlo3
          if(c==2) jlo=jlo-1
          if(c==3) klo=klo-1
          do k=klo,ixGhi3
            do j=jlo,ixGhi2
              do i=first,last,stride
                s(g)%ws(i,j,k,c)=s(g)%ws(anchor,j,k,c)
              end do
            end do
          end do
        end do
        do k=ixGlo3,ixGhi3
          do j=ixGlo2,ixGhi2
            do i=first,last,stride
              transverse=s(g)%surfaceC(i,j,k,2)*s(g)%ws(i,j,k,2) &
                -s(g)%surfaceC(i,j-1,k,2)*s(g)%ws(i,j-1,k,2) &
                +s(g)%surfaceC(i,j,k,3)*s(g)%ws(i,j,k,3) &
                -s(g)%surfaceC(i,j,k-1,3)*s(g)%ws(i,j,k-1,3)
              if(side==1) then
                s(g)%ws(i-1,j,k,1)=(s(g)%surfaceC(i,j,k,1)*s(g)%ws(i,j,k,1)+transverse) &
                  /s(g)%surfaceC(i-1,j,k,1)
              else
                s(g)%ws(i,j,k,1)=(s(g)%surfaceC(i-1,j,k,1)*s(g)%ws(i-1,j,k,1)-transverse) &
                  /s(g)%surfaceC(i,j,k,1)
              end if
            end do
          end do
        end do
      end do
    end do
    }
    {^IFONED call mpistop('CT radial fill requires 3D') }
    {^IFTWOD call mpistop('CT radial fill requires 3D') }
  end subroutine

  ! Caller supplies mhd_ct_sync_center. Every ghost cell is visited exactly once:
  ! x slabs, y slabs at active x, z slabs at active x/y. Angular slabs run first.
  ! Active total energy is
  ! unchanged; ghosts retain their prescribed thermal energy when B is replaced.
  subroutine yy_sync_ct_state(s,sync_center)
    type(state),intent(inout) :: s(max_blocks)
    interface
      subroutine sync_center(ixO^L,s,keep_pressure)
        use mod_global_parameters
        integer,intent(in) :: ixO^L
        type(state),intent(inout) :: s
        logical,intent(in) :: keep_pressure
      end subroutine
    end interface
    integer :: n,g,d,side,i,anchor,ixO^L,lo(ndim),hi(ndim),amin(ndim),amax(ndim),gmin(ndim),gmax(ndim)
    if(ndim/=3.or..not.yinyang_enabled.or.refine_max_level/=1.or.lb_automatic) &
      call mpistop('CT state sync requires fixed 3D Yin-Yang mesh')
    {^D& amin(^D)=ixMlo^D;amax(^D)=ixMhi^D;gmin(^D)=ixGlo^D;gmax(^D)=ixGhi^D;\}
    do n=1,igridstail
      g=igrids(n)
      call sync_center(ixM^LL,s(g),.false.)
      do d=ndim,1,-1
        do side=1,2
          lo=gmin;hi=gmax
          lo(1:d-1)=amin(1:d-1);hi(1:d-1)=amax(1:d-1)
          if(side==1) then
            hi(d)=amin(d)-1
          else
            lo(d)=amax(d)+1
          end if
          {^IFTHREED
          if(d==1) then
            anchor=0
            if(side==1.and.node(pig1_,g)==1) anchor=ixMlo1
            if(side==2.and.node(pig1_,g)==domain_nx1/block_nx1) anchor=ixMhi1
            if(anchor/=0) then
              ! Copy after active and angular ghost B/energy are consistent.
              do i=lo(1),hi(1)
                s(g)%w(i,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:)=s(g)%w(anchor,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:)
              end do
            end if
          end if
          }
          ixOmin^D=lo(^D);ixOmax^D=hi(^D);
          call sync_center(ixO^L,s(g),.true.)
        end do
      end do
    end do
  end subroutine

  ! Cache all owner angular-boundary edges, including radial endpoints.
  subroutine yy_runtime_ct_prepare(s)
    type(state),intent(in) :: s(max_blocks)
    {^IFTHREED
    integer,allocatable :: receivers(:,:)
    integer :: pass,n,g,c,i,j,k,idx(3),gi(3),bc(3),lo(3),hi(3),bn(3),nn(3),count
    if(runtime_patch%mesh_epoch==particle_mesh_epoch) return
    call yy_prepare_ct_patch(runtime_patch,s)
    call yy_prepare_ct_ghosts(runtime_ghosts,s,face_moments=yy_ct_face_moments)
    bn=[block_nx1,block_nx2,block_nx3];nn=[domain_nx1,domain_nx2,domain_nx3]
    hi=[ixMhi1,ixMhi2,ixMhi3]
    do pass=1,2
      count=0
      do n=1,igridstail
        g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
        do c=1,3
          lo=[ixMlo1,ixMlo2,ixMlo3]-1;lo(c)=lo(c)+1
          do k=lo(3),hi(3)
            do j=lo(2),hi(2)
              do i=lo(1),hi(1)
                idx=[i,j,k];gi=bc*bn+idx-[ixMlo1,ixMlo2,ixMlo3]
                if(.not.((c/=2.and.(gi(2)==-1.or.gi(2)==nn(2)-1)).or. &
                  (c/=3.and.(gi(3)==-1.or.gi(3)==nn(3)-1)))) cycle
                if(.not.yy_ct_edge_is_owner(g,idx,c)) cycle
                count=count+1
                if(pass==2) receivers(:,count)=[g,idx,c]
              end do
            end do
          end do
        end do
      end do
      if(pass==1) allocate(receivers(5,count))
    end do
    call yy_prepare_ct_edges(runtime_edges,s,receivers,radial_extrapolate=.true.)
    }
  end subroutine

  subroutine yy_runtime_ct_edges(s,fluxes)
    type(state),intent(in) :: s(max_blocks)
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    call yy_runtime_ct_prepare(s)
    call yy_ordered_ct_edges(runtime_patch,runtime_edges,fluxes)
  end subroutine

  subroutine yy_runtime_ct_ghosts(s)
    type(state),intent(inout) :: s(max_blocks)
    if(.not.associated(phys_ct_sync_center)) call mpistop('CT variable sync procedure is missing')
    call yy_runtime_ct_prepare(s)
    call yy_complete_ct_faces(runtime_patch,runtime_ghosts,s)
    call yy_sync_ct_state(s,phys_ct_sync_center)
  end subroutine

  ! Supply radial cell donors before angular primitive interpolation. Final
  ! radial B/energy is recomputed only after all face and angular data are ready.
  subroutine yy_ct_radial_cells(s)
    type(state),intent(inout) :: s(max_blocks)
    {^IFTHREED
    integer :: n,g,side,anchor,lo,hi,i
    do n=1,igridstail
      g=igrids(n)
      do side=1,2
        if(side==1.and.node(pig1_,g)/=1) cycle
        if(side==2.and.node(pig1_,g)/=domain_nx1/block_nx1) cycle
        anchor=ixMlo1;lo=ixGlo1;hi=ixMlo1-1
        if(side==2) then
          anchor=ixMhi1;lo=ixMhi1+1;hi=ixGhi1
        end if
        do i=lo,hi
          s(g)%w(i,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:)=s(g)%w(anchor,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:)
        end do
      end do
    end do
    }
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
    if(yy_ct_face_moments.and..not.yy_ct_experimental)call mpistop('face moments require experimental YY CT')
    if(yy_ct_experimental.and..not.yinyang_enabled) call mpistop('YY CT requires yinyang coordinates')
    if(.not.yinyang_enabled) return
    {^IFTHREED
    if(yy_ct_experimental) then
      if(.not.stagger_grid.or.physics_type/='mhd'.or..not.ct_stage_split.or. &
        .not.associated(phys_ct_sync_center).or..not.associated(phys_get_ct_edges).or. &
        .not.associated(phys_apply_ct_edges)) call mpistop('YY CT needs staged classical total-energy UCT-HLL')
      if(dimsplit.or.t_stepper/=threestep.or.t_integrator/=ssprk3.or.internalboundary) &
        call mpistop('YY CT needs unsplit SSPRK3; no internal boundaries')
      if(any(typeboundary(:,1:2)/=bc_cont)) call mpistop('YY CT requires continuous radial boundaries')
    else if(stagger_grid) then
      call mpistop('Yin-Yang CT requires explicit yy_ct_experimental=T')
    end if
    if(refine_max_level/=1.or.lb_automatic.or.local_timestep) &
      call mpistop('yinyang v1 requires fixed level-one mesh and global time stepping')
    if(B0field.or.use_particles.or.use_multigrid) &
      call mpistop('yinyang v1: B0, particles and multigrid are unsupported')
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
