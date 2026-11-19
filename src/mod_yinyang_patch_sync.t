! Cached same-patch copies on a fixed level-one grid. No averaging of EMFs.
module mod_yinyang_patch_sync
  use mpi
  implicit none
  private
  integer,parameter,public :: yy_sync_edges=1,yy_sync_faces=2,yy_sync_extended_faces=3
  public :: yy_patch_copy_plan,yy_patch_copy_build,yy_patch_copy_exchange
  type yy_patch_copy_plan
    integer :: comm=MPI_COMM_NULL
    integer,allocatable :: receivers(:,:),request(:,:),sc(:),rc(:),sd(:),rd(:)
    double precision,allocatable :: send(:),recv(:)
  end type
contains
  subroutine fail(comm,message)
    integer,intent(in) :: comm
    character(len=*),intent(in) :: message
    integer :: rank,ierr
    call MPI_COMM_RANK(comm,rank,ierr)
    write(*,'(a,i0,2a)') 'Yin-Yang patch sync rank ',rank,': ',message
    call MPI_ABORT(comm,51,ierr)
  end subroutine

  ! Directory: zero-based block coordinates, local ID, world rank, patch.
  ! Receivers/requests: local ID, zero-based local cell indices, component.
  ! A staggered upper face/edge at an internal interface belongs to the lower
  ! block in each transverse direction. Ownership never depends on MPI rank.
  subroutine yy_patch_copy_build(plan,comm,patch,kind,bn,nn,halo,directory)
    type(yy_patch_copy_plan),intent(out) :: plan
    integer,intent(in) :: comm,patch,kind,bn(3),nn(3),halo,directory(:,:)
    integer :: np,rank,ierr,nb(3),n,p,c,d,i,j,k,pass,nr,slot,idx(3),gi(3),bc(3),src(3),didx(3)
    integer :: lo(3),hi(3),amin(3),global_min(3),global_max(3),source_id,source_rank
    integer,allocatable :: owners(:,:,:,:),ids(:,:,:,:),cursor(:),out(:,:)
    logical :: faceaxis(3)
    call MPI_COMM_SIZE(comm,np,ierr);call MPI_COMM_RANK(comm,rank,ierr)
    if(kind<1.or.kind>3.or.patch<0.or.patch>1) call fail(comm,'invalid mode or patch')
    if(any(bn<=0).or.any(nn<=0).or.halo<0) call fail(comm,'invalid grid dimensions')
    if(any(mod(nn,bn)/=0).or.size(directory,1)/=6) call fail(comm,'invalid block directory shape')
    nb=nn/bn;plan%comm=comm
    allocate(owners(0:nb(1)-1,0:nb(2)-1,0:nb(3)-1,0:1),ids(0:nb(1)-1,0:nb(2)-1,0:nb(3)-1,0:1))
    owners=-1;ids=-1
    do n=1,size(directory,2)
      bc=directory(1:3,n);p=directory(6,n)
      if(any(bc<0).or.any(bc>=nb).or.p<0.or.p>1) call fail(comm,'invalid directory coordinate')
      if(directory(4,n)<1.or.directory(5,n)<0.or.directory(5,n)>=np) call fail(comm,'invalid block owner')
      if(owners(bc(1),bc(2),bc(3),p)/=-1) call fail(comm,'duplicate block coordinate')
      owners(bc(1),bc(2),bc(3),p)=directory(5,n);ids(bc(1),bc(2),bc(3),p)=directory(4,n)
    end do
    if(any(owners<0)) call fail(comm,'incomplete level-one mesh')
    allocate(plan%sc(np),plan%rc(np),plan%sd(np),plan%rd(np),cursor(np));plan%sc=0
    do pass=1,2
      nr=0
      do n=1,size(directory,2)
        if(directory(5,n)/=rank) cycle
        if(directory(6,n)/=patch) call fail(comm,'rank owns blocks in both patches')
        bc=directory(1:3,n)
        do c=1,3
          faceaxis=.false.;faceaxis(c)=.true.
          if(kind==yy_sync_edges) faceaxis=.not.faceaxis
          amin=0
          where(faceaxis) amin=-1
          lo=amin;hi=bn-1
          if(kind/=yy_sync_edges) then
            lo=lo-halo;hi=hi+halo
          end if
          global_min=amin;global_max=nn-1
          if(kind==yy_sync_extended_faces) then
            global_min(2:3)=global_min(2:3)-halo
            global_max(2:3)=global_max(2:3)+halo
          end if
          do k=lo(3),hi(3)
            do j=lo(2),hi(2)
              do i=lo(1),hi(1)
                idx=[i,j,k];gi=bc*bn+idx
                if(kind/=yy_sync_edges) then
                  ! Never copy active faces, even duplicated boundary anchors.
                  if(all(idx>=amin).and.all(idx<bn)) cycle
                end if
                if(any(gi<global_min).or.any(gi>global_max)) cycle
                src=max(0,min(nb-1,gi/bn))
                source_rank=owners(src(1),src(2),src(3),patch)
                source_id=ids(src(1),src(2),src(3),patch)
                if(source_rank==rank.and.source_id==directory(4,n)) cycle
                didx=gi-src*bn
                if(any(didx<lo).or.any(didx>hi)) call fail(comm,'source outside stored stencil')
                p=source_rank+1;nr=nr+1
                if(pass==1) then
                  plan%sc(p)=plan%sc(p)+1
                else
                  slot=plan%sd(p)+cursor(p)+1;cursor(p)=cursor(p)+1
                  out(:,slot)=[source_id,didx,c]
                  plan%receivers(:,slot)=[directory(4,n),idx,c]
                end if
              end do
            end do
          end do
        end do
      end do
      if(pass==1) then
        plan%sd(1)=0
        do p=2,np
          plan%sd(p)=plan%sd(p-1)+plan%sc(p-1)
        end do
        allocate(out(5,nr),plan%receivers(5,nr));cursor=0
      end if
    end do
    call MPI_ALLTOALL(plan%sc,1,MPI_INTEGER,plan%rc,1,MPI_INTEGER,comm,ierr)
    plan%rd(1)=0
    do p=2,np
      plan%rd(p)=plan%rd(p-1)+plan%rc(p-1)
    end do
    allocate(plan%request(5,sum(plan%rc)),plan%send(sum(plan%rc)),plan%recv(sum(plan%sc)))
    call MPI_ALLTOALLV(out,5*plan%sc,5*plan%sd,MPI_INTEGER,plan%request,5*plan%rc,5*plan%rd, &
      MPI_INTEGER,comm,ierr)
  end subroutine

  ! The caller packs every source before invoking this collective, then commits
  ! every received copy. Local-peer requests use the same ordering as MPI peers.
  subroutine yy_patch_copy_exchange(plan)
    type(yy_patch_copy_plan),intent(inout) :: plan
    integer :: ierr
    call MPI_ALLTOALLV(plan%send,plan%rc,plan%rd,MPI_DOUBLE_PRECISION, &
      plan%recv,plan%sc,plan%sd,MPI_DOUBLE_PRECISION,plan%comm,ierr)
  end subroutine
end module
