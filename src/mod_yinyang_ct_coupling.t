! CT coupling kernels. No implicit registration or change to CT enablement.
module mod_yinyang_ct_coupling
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mod_yinyang_staggered, only: yy_staggered_table, yy_staggered_build, &
    yy_staggered_exchange, yy_staggered_basis, yy_edge_location, yy_face_location
  implicit none
  private
  public :: yy_ct_edge_plan, yy_ct_edge_build, yy_ct_edge_exchange, yy_ct_edge_clear
  public :: yy_ct_edge_quadrature, yy_ct_close_ghosts
  public :: yy_ct_face_plan, yy_ct_face_build, yy_ct_face_exchange
  type yy_ct_face_plan
    type(yy_staggered_table) :: radial, tangent
    integer,allocatable :: rslot(:),tslot(:)
    double precision,allocatable :: normals(:,:),rvector(:,:),tvector(:,:)
    logical,allocatable :: covered(:)
  end type
  type yy_ct_edge_plan
    type(yy_staggered_table) :: samples
    double precision,allocatable :: line_weights(:,:,:), donor_length(:), scaled(:), vectors(:,:)
    logical,allocatable :: covered(:)
  end type
contains
  ! Br is patch-invariant as a scalar; tangential components require rotation.
  ! This split permits radial wall faces without borrowing physical ghost data.
  subroutine yy_ct_face_build(plan,comm,patch,bn,nn,radial,angular,directory,centers,directions,widths)
    type(yy_ct_face_plan),intent(out) :: plan
    integer,intent(in) :: comm,patch,bn(3),nn(3),directory(:,:),directions(:)
    double precision,intent(in) :: radial(0:),angular(2,2),centers(:,:)
    double precision,optional,intent(in) :: widths(:,:)
    double precision,allocatable :: rp(:,:),tp(:,:)
    double precision :: frame(3,3)
    integer :: n,nr,nt,ir,it,c
    if(size(centers,1)/=3.or.size(centers,2)/=size(directions)) &
      call edge_error(comm,'invalid face receiver shapes')
    if(any(directions<1).or.any(directions>3)) call edge_error(comm,'invalid face direction')
    nr=count(directions==1);nt=size(directions)-nr
    allocate(rp(3,nr),tp(3,nt),plan%rslot(nr),plan%tslot(nt))
    allocate(plan%normals(3,size(directions)),plan%covered(size(directions)))
    allocate(plan%rvector(3,nr),plan%tvector(3,nt))
    ir=0;it=0
    do n=1,size(directions)
      c=directions(n)
      call yy_staggered_basis(centers(2,n),centers(3,n),patch,frame)
      plan%normals(:,n)=frame(:,c)
      if(c==1) then
        ir=ir+1;rp(:,ir)=centers(:,n);plan%rslot(ir)=n
      else
        it=it+1;tp(:,it)=centers(:,n);plan%tslot(it)=n
      end if
    end do
    if(present(widths))then
      call yy_staggered_build(plan%radial,comm,patch,yy_face_location,bn,nn,radial,angular,directory,rp, &
        face_widths=widths(:,plan%rslot),face_directions=directions(plan%rslot))
      call yy_staggered_build(plan%tangent,comm,patch,yy_face_location,bn,nn,radial,angular,directory,tp, &
        face_widths=widths(:,plan%tslot),face_directions=directions(plan%tslot))
    else
      call yy_staggered_build(plan%radial,comm,patch,yy_face_location,bn,nn,radial,angular,directory,rp, &
        radial_scalar=.true.)
      call yy_staggered_build(plan%tangent,comm,patch,yy_face_location,bn,nn,radial,angular,directory,tp)
    end if
    plan%covered(plan%rslot)=plan%radial%covered
    plan%covered(plan%tslot)=plan%tangent%covered
  end subroutine

  subroutine yy_ct_face_exchange(plan,radial_samples,tangent_samples,values)
    type(yy_ct_face_plan),intent(inout) :: plan
    double precision,intent(in) :: radial_samples(:),tangent_samples(:)
    double precision,intent(out) :: values(:)
    integer :: n,t
    if(size(values)/=size(plan%covered)) call edge_error(plan%radial%comm,'invalid face output size')
    call yy_staggered_exchange(plan%radial,radial_samples,plan%rvector)
    call yy_staggered_exchange(plan%tangent,tangent_samples,plan%tvector)
    values=0.d0
    do n=1,size(plan%rslot)
      t=plan%rslot(n)
      if(plan%covered(t)) values(t)=dot_product(plan%normals(:,t),plan%rvector(:,n))
    end do
    do n=1,size(plan%tslot)
      t=plan%tslot(n)
      if(plan%covered(t)) values(t)=dot_product(plan%normals(:,t),plan%tvector(:,n))
    end do
  end subroutine

  subroutine edge_error(comm,message)
    integer,intent(in) :: comm
    character(len=*),intent(in) :: message
    integer :: rank,ierr
    call MPI_COMM_RANK(comm,rank,ierr)
    write(*,'(a,i0,2a)') 'Yin-Yang CT coupling rank ',rank,': ',message
    call MPI_ABORT(comm,31,ierr)
  end subroutine

  subroutine yy_ct_edge_clear(plan)
    type(yy_ct_edge_plan),intent(inout) :: plan
    type(yy_ct_edge_plan) :: empty
    plan=empty
  end subroutine

  ! Positive-coordinate oriented edge, centered at position. Width is the
  ! coordinate width along direction (dr, dtheta or dphi), not arc length.
  pure subroutine yy_ct_edge_quadrature(patch,direction,position,width,points,weights)
    integer,intent(in) :: patch,direction
    double precision,intent(in) :: position(3),width
    double precision,intent(out) :: points(3,2),weights(3,2)
    double precision :: frame(3,3),jacobian
    integer :: q
    do q=1,2
      points(:,q)=position
      points(direction,q)=position(direction)+(2*q-3)*width/(2.d0*sqrt(3.d0))
      call yy_staggered_basis(points(2,q),points(3,q),patch,frame)
      jacobian=1.d0
      if(direction==2) jacobian=points(1,q)
      if(direction==3) jacobian=points(1,q)*sin(points(2,q))
      weights(:,q)=0.5d0*width*jacobian*frame(:,direction)
    end do
  end subroutine

  subroutine yy_ct_edge_build(plan,comm,patch,bn,nn,radial,angular,directory,centers,widths,directions,radial_extrapolate,native_midpoint)
    type(yy_ct_edge_plan),intent(inout) :: plan
    integer,intent(in) :: comm,patch,bn(3),nn(3),directory(:,:),directions(:)
    double precision,intent(in) :: radial(0:),angular(2,2),centers(:,:),widths(:)
    logical,optional,intent(in) :: radial_extrapolate,native_midpoint
    double precision,allocatable :: points(:,:)
    integer,allocatable :: idmap(:)
    double precision :: pos(3,2),weights(3,2),r,step(2),frame(3,3),jacobian
    integer :: n,nt,c,rank,ierr,g,ig,idx,ns,nq
    logical :: midpoint
    call yy_ct_edge_clear(plan)
    nt=size(directions)
    if(size(centers,1)/=3.or.size(centers,2)/=nt.or.size(widths)/=nt) &
      call edge_error(comm,'invalid receiver array shapes')
    if(any(directions<1).or.any(directions>3).or.any(widths<=0.d0).or. &
      .not.all(ieee_is_finite(widths))) call edge_error(comm,'invalid edge direction or width')
    ! Native UCT supplies midpoint approximations, not exact line averages.
    ! Select quadrature once at plan construction; the exchange kernel is shared.
    midpoint=.false.;if(present(native_midpoint))midpoint=native_midpoint
    nq=2;if(midpoint)nq=1
    allocate(points(3,nq*nt),plan%line_weights(3,nq,nt),plan%covered(nt))
    do n=1,nt
      if(midpoint)then
        c=directions(n)
        call yy_staggered_basis(centers(2,n),centers(3,n),patch,frame)
        jacobian=1.d0
        if(c==2)jacobian=centers(1,n)
        if(c==3)jacobian=centers(1,n)*sin(centers(2,n))
        points(:,n)=centers(:,n)
        plan%line_weights(:,1,n)=widths(n)*jacobian*frame(:,c)
      else
        call yy_ct_edge_quadrature(patch,directions(n),centers(:,n),widths(n),pos,weights)
        points(:,2*n-1:2*n)=pos;plan%line_weights(:,:,n)=weights
      end if
    end do
    call yy_staggered_build(plan%samples,comm,patch,yy_edge_location,bn,nn,radial,angular,directory,points, &
      edge_radial_extrapolate=radial_extrapolate)
    ns=size(plan%samples%request,2)
    allocate(plan%donor_length(ns),plan%scaled(ns),plan%vectors(3,nq*nt))
    call MPI_COMM_RANK(comm,rank,ierr)
    allocate(idmap(maxval(directory(4,:))));idmap=0
    do n=1,size(directory,2)
      if(directory(5,n)==rank) idmap(directory(4,n))=n
    end do
    step=(angular(:,2)-angular(:,1))/nn(2:3)
    do n=1,ns
      g=plan%samples%request(1,n);ig=idmap(g);c=plan%samples%request(5,n)
      idx=directory(1,ig)*bn(1)+plan%samples%request(2,n)
      if(c==1) then
        plan%donor_length(n)=radial(idx+1)-radial(idx)
      else
        r=radial(idx+1)
        if(c==2) then
          plan%donor_length(n)=r*step(1)
        else
          ! xy magnitude is invariant under Yin-Yang rotation only after
          ! converting back to the donor patch; use its angular index instead.
          idx=directory(2,ig)*bn(2)+plan%samples%request(3,n)
          plan%donor_length(n)=r*sin(angular(1,1)+(idx+1)*step(1))*step(2)
        end if
      end if
    end do
    if(any(plan%donor_length<=0.d0)) call edge_error(comm,'nonpositive donor edge length')
    do n=1,nt
      plan%covered(n)=all(plan%samples%covered(nq*(n-1)+1:nq*n))
    end do
  end subroutine

  ! Donor input is the native dt*length*E_tangent. Divide only by length,
  ! reconstruct dt*E, then apply the plan quadrature. Runtime UCT uses midpoint
  ! quadrature on BOTH donors and receivers. Zero dt is safe.
  ! All ranks call, including those with zero receivers. No runtime allocation.
  subroutine yy_ct_edge_exchange(plan,donor_integrals,receiver_integrals)
    type(yy_ct_edge_plan),intent(inout) :: plan
    double precision,intent(in) :: donor_integrals(:)
    double precision,intent(out) :: receiver_integrals(:)
    integer :: n,q,nq
    nq=size(plan%line_weights,2)
    if(size(donor_integrals)/=size(plan%scaled).or.size(receiver_integrals)/=size(plan%covered)) &
      call edge_error(plan%samples%comm,'edge exchange shapes disagree with plan')
    plan%scaled=donor_integrals/plan%donor_length
    call yy_staggered_exchange(plan%samples,plan%scaled,plan%vectors)
    receiver_integrals=0.d0
    do n=1,size(plan%covered)
      if(.not.plan%covered(n)) cycle
      do q=1,nq
        receiver_integrals(n)=receiver_integrals(n)+ &
          dot_product(plan%line_weights(:,q,n),plan%vectors(:,nq*(n-1)+q))
      end do
    end do
  end subroutine

  ! Face-index convention: face(i,j,k,c) is the upper c-face of cell(i,j,k).
  ! Only angular ghost cells at selected physical patch sides are constrained.
  ! theta slabs first (phi interior), then phi slabs including theta corners.
  ! Existing tangential ghost faces and active-boundary normal faces are inputs.
  ! Radial physical ghosts and same-patch inter-block ghosts are not modified.
  subroutine yy_ct_close_ghosts(array_lo,active_lo,active_hi,ghost_lo,ghost_hi,sides,faces,area,status)
    integer,intent(in) :: array_lo(3),active_lo(3),active_hi(3),ghost_lo(3),ghost_hi(3)
    logical,intent(in) :: sides(2,2)
    double precision,intent(inout) :: faces(array_lo(1):,array_lo(2):,array_lo(3):,:)
    double precision,intent(in) :: area(array_lo(1):,array_lo(2):,array_lo(3):,:)
    integer,intent(out) :: status
    integer :: i,j,k,side,jlo,jhi,start,finish,stride,lo(3),hi(3),upper(4),pass
    double precision :: transverse,flux
    status=1
    lo=ghost_lo;hi=ghost_hi
    if(any(active_lo>active_hi).or.any(lo>active_lo).or.any(hi<active_hi)) return
    if(any(array_lo>lo-1).or.any(shape(faces)/=shape(area))) return
    upper=ubound(faces)
    if(any(hi>upper(1:3))) return
    if(size(faces,4)/=3) return
    ! Validate the complete input footprint before any writes (atomic failure).
    do pass=1,2
    do side=1,2
      if(.not.sides(side,1)) cycle
      start=active_lo(2)-1;finish=lo(2);stride=-1
      if(side==2) then
        start=active_hi(2)+1;finish=hi(2);stride=1
      end if
      do j=start,finish,stride
        do k=active_lo(3),active_hi(3)
          do i=active_lo(1),active_hi(1)
            if(pass==1) then
              if(.not.valid_cell(i,j,k)) return
              cycle
            end if
            transverse=area(i,j,k,1)*faces(i,j,k,1)-area(i-1,j,k,1)*faces(i-1,j,k,1) &
              +area(i,j,k,3)*faces(i,j,k,3)-area(i,j,k-1,3)*faces(i,j,k-1,3)
            if(side==1) then
              flux=area(i,j,k,2)*faces(i,j,k,2)+transverse
              faces(i,j-1,k,2)=flux/area(i,j-1,k,2)
            else
              flux=area(i,j-1,k,2)*faces(i,j-1,k,2)-transverse
              faces(i,j,k,2)=flux/area(i,j,k,2)
            end if
          end do
        end do
      end do
    end do
    jlo=active_lo(2);jhi=active_hi(2)
    if(sides(1,1)) jlo=lo(2)
    if(sides(2,1)) jhi=hi(2)
    do side=1,2
      if(.not.sides(side,2)) cycle
      start=active_lo(3)-1;finish=lo(3);stride=-1
      if(side==2) then
        start=active_hi(3)+1;finish=hi(3);stride=1
      end if
      do k=start,finish,stride
        do j=jlo,jhi
          do i=active_lo(1),active_hi(1)
            if(pass==1) then
              if(.not.valid_cell(i,j,k)) return
              cycle
            end if
            transverse=area(i,j,k,1)*faces(i,j,k,1)-area(i-1,j,k,1)*faces(i-1,j,k,1) &
              +area(i,j,k,2)*faces(i,j,k,2)-area(i,j-1,k,2)*faces(i,j-1,k,2)
            if(side==1) then
              flux=area(i,j,k,3)*faces(i,j,k,3)+transverse
              faces(i,j,k-1,3)=flux/area(i,j,k-1,3)
            else
              flux=area(i,j,k-1,3)*faces(i,j,k-1,3)-transverse
              faces(i,j,k,3)=flux/area(i,j,k,3)
            end if
          end do
        end do
      end do
    end do
    end do ! validation pass, then correction pass
    status=0
  contains
    logical function valid_cell(i,j,k) result(valid)
      integer,intent(in) :: i,j,k
      double precision :: areas(6),values(6)
      areas=[area(i,j,k,1),area(i-1,j,k,1),area(i,j,k,2), &
        area(i,j-1,k,2),area(i,j,k,3),area(i,j,k-1,3)]
      values=[faces(i,j,k,1),faces(i-1,j,k,1),faces(i,j,k,2), &
        faces(i,j-1,k,2),faces(i,j,k,3),faces(i,j,k-1,3)]
      valid=all(ieee_is_finite(areas)).and.all(ieee_is_finite(values)).and.all(areas>0.d0)
    end function
  end subroutine
end module
