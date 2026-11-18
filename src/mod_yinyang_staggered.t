! Geometry-only, fixed-grid staggered sampling. No solver state is updated.
module mod_yinyang_staggered
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: yy_face_location=1, yy_edge_location=2
  public :: yy_staggered_table, yy_staggered_build, yy_staggered_exchange
  public :: yy_staggered_clear, yy_staggered_basis, yy_staggered_offset
  type yy_staggered_table
    integer :: comm=MPI_COMM_NULL, stencil=24
    integer, allocatable :: request(:,:), slot(:,:), sc(:),rc(:),sd(:),rd(:)
    logical, allocatable :: covered(:)
    double precision, allocatable :: donor_position(:,:), donor_basis(:,:)
    double precision, allocatable :: coefficient(:,:,:), send(:),recv(:)
  end type
contains
  subroutine fail(comm,code)
    integer,intent(in) :: comm,code
    integer :: ierr,rank
    character(len=96) :: message
    select case(code)
    case(10);message='invalid patch or sampling family'
    case(11,12);message='invalid grid sizes or array shapes'
    case(13);message='invalid radial faces or angular extents'
    case(14,15);message='invalid block directory entry'
    case(16,17);message='duplicate or missing level-one block'
    case(18);message='singular or ill-conditioned Cartesian reconstruction'
    case(19,20);message='donor request routed to the wrong owner or patch'
    case(21);message='exchange arrays disagree with table'
    case(23);message='duplicate local block ID on one rank'
    case(24);message='empty global block directory'
    case(25);message='nonfinite coordinate input'
    case default;message='invalid staggered table'
    end select
    call MPI_COMM_RANK(comm,rank,ierr)
    write(*,'(a,i0,2a)') 'Yin-Yang staggered table rank ',rank,': ',trim(message)
    call MPI_ABORT(comm,code,ierr)
  end subroutine

  pure subroutine yy_staggered_offset(family,c,offset)
    integer,intent(in) :: family,c
    double precision,intent(out) :: offset(3)
    offset=0.5d0
    if(family==yy_face_location) then
      offset(c)=1.d0
    else
      offset=1.d0;offset(c)=0.5d0
    end if
  end subroutine

  pure subroutine yy_staggered_basis(theta,phi,patch,b)
    double precision,intent(in) :: theta,phi
    integer,intent(in) :: patch
    double precision,intent(out) :: b(3,3)
    double precision :: tmp(3)
    b(:,1)=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
    b(:,2)=[cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)]
    b(:,3)=[-sin(phi),cos(phi),0.d0]
    if(patch==1) then
      b(1,:)=-b(1,:);tmp=b(2,:);b(2,:)=b(3,:);b(3,:)=tmp
    end if
  end subroutine

  subroutine yy_staggered_clear(t)
    type(yy_staggered_table),intent(inout) :: t
    type(yy_staggered_table) :: empty
    t=empty
  end subroutine

  ! Directory rows: zero-based block coordinates, local block ID, world rank,
  ! patch (0/1). Every rank supplies the same complete level-one directory.
  ! Receiver points are (r,theta,phi) in its own patch. Donors are scalar normal
  ! face components or tangent edge components, NOT dt*length edge integrals.
  ! Request rows: local block ID, three zero-based local cell indices, component.
  ! An index -1 denotes the physical lower face/edge, not an unfilled ghost donor.
  subroutine yy_staggered_build(t,comm,patch,family,bn,nn,radial,angular,directory,points,radial_scalar)
    type(yy_staggered_table),intent(inout) :: t
    integer,intent(in) :: comm,patch,family,bn(3),nn(3),directory(:,:)
    double precision,intent(in) :: radial(0:),angular(2,2),points(:,:)
    logical,optional,intent(in) :: radial_scalar
    integer :: np,rank,ierr,n,p,c,q,a,b,d,h,pass,target,base(3),idx(3),bc(3),loc(3),components
    integer :: qrank(24),qreq(5,24),owner,nt,ns,nr,off,ib,nb(3),lo,hi,mid
    integer,allocatable :: owners(:,:,:,:),ids(:,:,:,:),cursor(:),out(:,:),idmap(:)
    double precision :: xyz(3),v(3),sph(3),step(2),offset(3),f(3),u,coord(3),frame(3,3)
    double precision :: w(24),basis(3,24),normal(3,3),inv(3,3),det,left,right,x,tolerance
    logical :: valid,scalar
    call yy_staggered_clear(t)
    scalar=.false.
    if(present(radial_scalar)) scalar=radial_scalar
    components=3
    if(scalar) then
      if(family/=yy_face_location) call fail(comm,10)
      components=1;t%stencil=8
    end if
    call MPI_COMM_SIZE(comm,np,ierr);call MPI_COMM_RANK(comm,rank,ierr)
    if(patch<0.or.patch>1.or.(family/=yy_face_location.and.family/=yy_edge_location)) &
      call fail(comm,10)
    if(any(bn<=0).or.any(nn<2)) call fail(comm,11)
    if(any(mod(nn,bn)/=0).or.size(radial)/=nn(1)+1.or.size(points,1)/=3.or.size(directory,1)/=6) &
      call fail(comm,12)
    if(.not.all(ieee_is_finite(radial)).or..not.all(ieee_is_finite(angular)).or. &
      .not.all(ieee_is_finite(points))) call fail(comm,25)
    if(any(radial(1:)<=radial(:nn(1)-1)).or.radial(0)<=0.d0.or.any(angular(:,2)<=angular(:,1))) &
      call fail(comm,13)
    if(size(directory,2)==0) call fail(comm,24)
    if(any(directory(4,:)<1)) call fail(comm,15)
    nb=nn/bn;t%comm=comm;nt=size(points,2)
    allocate(owners(0:nb(1)-1,0:nb(2)-1,0:nb(3)-1,0:1),ids(0:nb(1)-1,0:nb(2)-1,0:nb(3)-1,0:1))
    allocate(idmap(maxval(directory(4,:))))
    idmap=0;owners=-1;ids=-1
    do n=1,size(directory,2)
      bc=directory(1:3,n);p=directory(6,n)
      if(any(bc<0).or.any(bc>=nb).or.p<0.or.p>1) call fail(comm,14)
      if(directory(5,n)<0.or.directory(5,n)>=np.or.directory(4,n)<1) call fail(comm,15)
      if(owners(bc(1),bc(2),bc(3),p)/=-1) call fail(comm,16)
      if(directory(5,n)==rank) then
        if(idmap(directory(4,n))/=0) call fail(comm,23)
        idmap(directory(4,n))=n
      end if
      owners(bc(1),bc(2),bc(3),p)=directory(5,n);ids(bc(1),bc(2),bc(3),p)=directory(4,n)
    end do
    if(any(owners<0)) call fail(comm,17)
    allocate(t%sc(np),t%rc(np),t%sd(np),t%rd(np),cursor(np),t%covered(nt),t%slot(t%stencil,nt),t%coefficient(3,t%stencil,nt))
    t%sc=0;t%slot=0;t%coefficient=0.d0;t%covered=.false.
    step=(angular(:,2)-angular(:,1))/nn(2:3)
    do pass=1,2
      do target=1,nt
        call yy_staggered_basis(points(2,target),points(3,target),patch,frame)
        xyz=points(1,target)*frame(:,1);v=xyz
        if(patch==0) v=[-xyz(1),xyz(3),xyz(2)]
        sph=[points(1,target),atan2(sqrt(v(1)**2+v(2)**2),v(3)),atan2(v(2),v(1))]
        valid=.true.;q=0
        do c=1,components
          call yy_staggered_offset(family,c,offset)
          ! Bracket the shared nonuniform radial mesh without extrapolation.
          lo=0;hi=nn(1)-1
          if(offset(1)==1.d0) lo=-1
          left=radial_sample(lo,offset(1));right=radial_sample(hi,offset(1))
          x=sph(1)
          tolerance=1.d-12*max(1.d0,abs(left),abs(right))
          if(x<left-tolerance.or.x>right+tolerance) then
            valid=.false.;exit
          end if
          x=max(left,min(right,x))
          do while(hi-lo>1)
            mid=(lo+hi)/2
            if(radial_sample(mid,offset(1))<=x) then
              lo=mid
            else
              hi=mid
            end if
          end do
          base(1)=lo;left=radial_sample(lo,offset(1));right=radial_sample(lo+1,offset(1))
          f(1)=(x-left)/(right-left)
          do d=2,3
            u=(sph(d)-angular(d-1,1))/step(d-1)-offset(d)
            lo=0;if(offset(d)==1.d0) lo=-1
            if(u<lo-1.d-12.or.u>nn(d)-1+1.d-12) then
              valid=.false.;exit
            end if
            u=max(dble(lo),min(dble(nn(d)-1),u))
            base(d)=min(floor(u),nn(d)-2);f(d)=u-base(d)
          end do
          if(.not.valid) exit
          do h=0,1
            do b=0,1
              do a=0,1
                q=q+1;idx=base+[a,b,h];bc=max(0,idx/bn);loc=idx-bc*bn
                owner=owners(bc(1),bc(2),bc(3),1-patch)
                qrank(q)=owner+1;qreq(:,q)=[ids(bc(1),bc(2),bc(3),1-patch),loc,c]
                w(q)=merge(1.d0-f(1),f(1),a==0)*merge(1.d0-f(2),f(2),b==0)*merge(1.d0-f(3),f(3),h==0)
                coord=[radial_sample(idx(1),offset(1)),angular(:,1)+(idx(2:3)+offset(2:3))*step]
                call yy_staggered_basis(coord(2),coord(3),1-patch,frame)
                basis(:,q)=frame(:,c)
              end do
            end do
          end do
        end do
        if(.not.valid) cycle
        t%covered(target)=.true.
        if(pass==1) then
          do q=1,t%stencil
            t%sc(qrank(q))=t%sc(qrank(q))+1
          end do
        else
          ! Weighted Cartesian constant fit. Normalization makes any uniform
          ! Cartesian vector exact despite non-collocated component samples.
          if(.not.scalar) then
          normal=0.d0
          do q=1,t%stencil
            do d=1,3
              normal(:,d)=normal(:,d)+w(q)*basis(:,q)*basis(d,q)
            end do
          end do
          inv(:,1)=cross(normal(:,2),normal(:,3))
          inv(:,2)=cross(normal(:,3),normal(:,1))
          inv(:,3)=cross(normal(:,1),normal(:,2))
          det=dot_product(normal(:,1),inv(:,1))
          if(det<1.d-8) call fail(comm,18)
          inv=transpose(inv)/det
          end if
          do q=1,t%stencil
            p=qrank(q);off=t%sd(p)+cursor(p)+1;cursor(p)=cursor(p)+1
            out(:,off)=qreq(:,q);t%slot(q,target)=off
            if(scalar) then
              ! Br is a scalar under patch rotation at a fixed physical point.
              ! Its own face lattice reaches radial walls without extrapolation.
              t%coefficient(:,q,target)=w(q)*xyz/points(1,target)
            else
              t%coefficient(:,q,target)=w(q)*matmul(inv,basis(:,q))
            end if
          end do
        end if
      end do
      if(pass==1) then
        t%sd(1)=0
        do p=2,np
          t%sd(p)=t%sd(p-1)+t%sc(p-1)
        end do
        allocate(out(5,sum(t%sc)));cursor=0
      end if
    end do
    call MPI_ALLTOALL(t%sc,1,MPI_INTEGER,t%rc,1,MPI_INTEGER,comm,ierr)
    t%rd(1)=0
    do p=2,np
      t%rd(p)=t%rd(p-1)+t%rc(p-1)
    end do
    ns=sum(t%rc);nr=sum(t%sc)
    allocate(t%request(5,ns),t%donor_position(3,ns),t%donor_basis(3,ns),t%send(ns),t%recv(nr))
    call MPI_ALLTOALLV(out,5*t%sc,5*t%sd,MPI_INTEGER,t%request,5*t%rc,5*t%rd,MPI_INTEGER,comm,ierr)
    do n=1,ns
      ib=t%request(1,n);owner=idmap(ib)
      if(owner==0) call fail(comm,19)
      if(directory(6,owner)/=patch) call fail(comm,20)
      c=t%request(5,n);idx=directory(1:3,owner)*bn+t%request(2:4,n)
      call yy_staggered_offset(family,c,offset)
      coord=[radial_sample(idx(1),offset(1)),angular(:,1)+(idx(2:3)+offset(2:3))*step]
      call yy_staggered_basis(coord(2),coord(3),patch,frame)
      t%donor_position(:,n)=coord(1)*frame(:,1);t%donor_basis(:,n)=frame(:,c)
    end do
  contains
    pure double precision function radial_sample(i,fraction) result(r)
      integer,intent(in) :: i
      double precision,intent(in) :: fraction
      if(fraction==1.d0) then
        r=radial(i+1)
      else
        r=0.5d0*(radial(i)+radial(i+1))
      end if
    end function
  end subroutine

  ! Call once per rank outside block loops. Input samples follow request order;
  ! output is Cartesian. Uncovered receivers are zero and MUST be masked.
  subroutine yy_staggered_exchange(t,samples,cartesian)
    type(yy_staggered_table),intent(inout) :: t
    double precision,intent(in) :: samples(:)
    double precision,intent(out) :: cartesian(:,:)
    integer :: ierr,n,q
    if(size(samples)/=size(t%send).or.size(cartesian,1)/=3.or.size(cartesian,2)/=size(t%covered)) &
      call fail(t%comm,21)
    t%send=samples
    call MPI_ALLTOALLV(t%send,t%rc,t%rd,MPI_DOUBLE_PRECISION,t%recv,t%sc,t%sd,MPI_DOUBLE_PRECISION,t%comm,ierr)
    cartesian=0.d0
    do n=1,size(t%covered)
      if(.not.t%covered(n)) cycle
      do q=1,t%stencil
        cartesian(:,n)=cartesian(:,n)+t%coefficient(:,q,n)*t%recv(t%slot(q,n))
      end do
    end do
  end subroutine

  pure function cross(a,b) result(c)
    double precision,intent(in) :: a(3),b(3)
    double precision :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function
end module
