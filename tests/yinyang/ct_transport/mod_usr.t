! Full MHD weak compact-loop advection through experimental Yin-Yang CT.
module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved,phys_get_pthermal
  use mod_yinyang, only: yy_basis,yy_position
  use mod_functions_bfield, only: get_divb
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical :: initialized=.false.
  logical :: native_control=.false.,save_samples=.false.
  integer :: reference_patch=0,frame_patch=0
  integer :: direction=1,placement=1
  double precision :: amplitude=1.d-3,radius=0.35d0,origin(3),velocity(3),normal(3),a0
  integer :: stage_calls=0,last_iteration=-1,history_unit,test_world_rank
  double precision :: largest_emf=0.d0
contains
  subroutine usr_init()
    integer :: n
    double precision :: phi
    namelist /transport_test/ direction,placement,amplitude,radius,native_control,reference_patch,save_samples
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,transport_test,end=10)
10    close(unitpar)
    end do
    call MPI_COMM_RANK(MPI_COMM_WORLD,test_world_rank,ierrmpi)
    if(native_control)then
      call set_coordinate_system('spherical')
    else
      call set_coordinate_system('yinyang')
    end if
    phi=0.d0
    if(placement==1)phi=146.25d0*acos(-1.d0)/180.d0
    velocity=direction*[-sin(phi),cos(phi),0.d0]
    normal=[cos(phi),sin(phi),0.d0]
    origin=1.5d0*normal-0.3d0*velocity
    a0=amplitude*radius/((8.d0/sqrt(7.d0))*(6.d0/7.d0)**3)
    if(radius<=0.d0.or.amplitude<=0.d0)call mpistop('invalid loop parameters')
    call mhd_activate()
    if(native_control)usr_special_bc=>native_fluid_boundary
    usr_init_one_grid=>initial
    usr_init_vector_potential=>potential
    usr_ct_stage=>stage_diagnostic
    usr_process_global=>diagnostic
    usr_process_adv_global=>advanced
  end subroutine
  ! Prescribe only the uniform fluid in angular physical ghosts; B uses cont.
  ! This avoids imposing constant spherical components on a Cartesian flow.
  subroutine native_fluid_boundary(qdt,qt,ixI^L,ixO^L,iB,w,x)
    double precision,intent(in) :: qdt,qt,x(ixI^S,ndim)
    integer,intent(in) :: ixI^L,ixO^L,iB
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: ix^D
    double precision :: frame(3,3)
    if(iB<=2)call mpistop('native fluid boundary is angular only')
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=0.1d0
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call yy_basis(x(ix^D,2),x(ix^D,3),reference_patch,frame)
      w(ix^D,mom(:))=matmul(transpose(frame),velocity)
    {end do\}
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: ix^D
    double precision :: frame(3,3)
    frame_patch=yy_patch
    if(native_control)frame_patch=reference_patch
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=0.1d0
    frame_patch=yy_patch
    if(native_control)frame_patch=reference_patch
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call yy_basis(x(ix^D,2),x(ix^D,3),frame_patch,frame)
      w(ix^D,mom(:))=matmul(transpose(frame),velocity)
    {end do\}
    call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
    call mhd_face_to_center(ixO^L,block)
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine
  subroutine potential(ixI^L,ixO^L,x,A,idir)
    integer,intent(in) :: ixI^L,ixO^L,idir
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(out) :: A(ixI^S)
    integer :: ix^D
    double precision :: frame(3,3),pos(3),cart(3),q(3),ss
    frame_patch=yy_patch
    if(native_control)frame_patch=reference_patch
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call yy_basis(x(ix^D,2),x(ix^D,3),frame_patch,frame)
      call yy_position(x(ix^D,1),x(ix^D,2),x(ix^D,3),frame_patch,pos)
      q=pos-origin
      ss=sum(q*q)/radius**2
      cart=0.d0
      if(ss<1.d0)cart=a0*(1.d0-ss)**4*normal
      A(ix^D)=dot_product(cart,frame(:,idir))
    {end do\}
  end subroutine
  subroutine stage_diagnostic(qt,qdt,sCT,s,fluxes)
    double precision,intent(in) :: qt,qdt
    type(state),intent(in) :: sCT(max_blocks),s(max_blocks)
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    integer :: n,g
    stage_calls=stage_calls+1
    do n=1,igridstail
      g=igrids(n)
      if(.not.all(ieee_is_finite(fluxes(g)%edge))) call mpistop('nonfinite static-test EMF')
      largest_emf=max(largest_emf,maxval(abs(fluxes(g)%edge)))
    end do
  end subroutine
  subroutine advanced(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    call diagnostic(iit+1,qt+dt)
  end subroutine
  subroutine diagnostic(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    integer :: n,g,i,j,k,c,bc(3),gi(3),bn(3),nn(3),lo(3),hi(3),fl(3),u,callsmin,callsmax
    double precision :: divb(ixG^T),pth(ixG^T),local(11),maximum(11),mins(2),minimum(2),tot(12),total(12)
    double precision :: error,frame(3,3),expected(3),pos(3),cart(3),q(3),ss,weight,other(3),theta,phi,vol,b2
    logical :: seam,radial,active
    character(len=80) :: name
    if(iit==last_iteration) return
    last_iteration=iit
    frame_patch=yy_patch
    if(native_control)frame_patch=reference_patch
    if(.not.initialized) then
      if(iit/=0)call mpistop('transport test missing initial record')
      initialized=.true.
      if(test_world_rank==0)open(newunit=history_unit,file='transport-history.txt',status='replace')
      call sleep(1)
    end if
    local=0.d0;mins=huge(1.d0);tot=0.d0
    lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3]
    bn=[block_nx1,block_nx2,block_nx3];nn=[domain_nx1,domain_nx2,domain_nx3]
    do n=1,igridstail
      g=igrids(n);block=>ps(g)
      bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
      if(.not.all(ieee_is_finite(ps(g)%w)))call mpistop('nonfinite static-test state')
      call get_divb(ps(g)%w,ixG^LL,ixG^LL,divb)
      call phys_get_pthermal(ps(g)%w,ps(g)%x,ixG^LL,ixM^LL,pth)
      do k=ixGlo3,ixGhi3
        do j=ixGlo2,ixGhi2
          do i=ixGlo1,ixGhi1
            gi=bc*bn+[i,j,k]-lo
            active=all([i,j,k]>=lo).and.all([i,j,k]<=hi)
            radial=gi(1)<0.or.gi(1)>=nn(1)
            if(radial)local(4)=max(local(4),abs(divb(i,j,k)))
            if(.not.active)then
              local(3)=max(local(3),abs(divb(i,j,k)))
              cycle
            end if
            seam=any(gi(2:3)<2).or.any(gi(2:3)>=nn(2:3)-2)
            local(1)=max(local(1),abs(divb(i,j,k)))
            if(seam)local(2)=max(local(2),abs(divb(i,j,k)))
            call yy_basis(ps(g)%x(i,j,k,2),ps(g)%x(i,j,k,3),frame_patch,frame)
            call yy_position(ps(g)%x(i,j,k,1),ps(g)%x(i,j,k,2),ps(g)%x(i,j,k,3),frame_patch,pos)
            q=pos-origin-velocity*qt
            ss=sum(q*q)/radius**2
            expected=0.d0
            if(ss<1.d0)expected=-8.d0*a0/radius**2*(1.d0-ss)**3* &
              [q(2)*normal(3)-q(3)*normal(2),q(3)*normal(1)-q(1)*normal(3),q(1)*normal(2)-q(2)*normal(1)]
            cart=matmul(frame,ps(g)%w(i,j,k,mag(:)))
            error=sqrt(sum((cart-expected)**2))
            local(5)=max(local(5),error)
            local(6)=max(local(6),sqrt(sum((matmul(frame,ps(g)%w(i,j,k,mom(:)))/ps(g)%w(i,j,k,rho_)-velocity)**2)))
            if(seam)local(7)=max(local(7),error)
            if(.not.seam.and.gi(1)>=2.and.gi(1)<nn(1)-2)local(8)=max(local(8),error)
            local(10)=max(local(10),sqrt(sum(cart*cart)))
            local(11)=max(local(11),abs(ps(g)%w(i,j,k,rho_)-1.d0))
            ! Partition of unity: weight 1/2 where both patches cover this point.
            other=pos
            if(frame_patch==0)other=[-pos(1),pos(3),pos(2)]
            theta=acos(other(3)/sqrt(sum(other*other)));phi=atan2(other(2),other(1))
            weight=1.d0
            if(.not.native_control)then
              if(theta>=xprobmin2.and.theta<=xprobmax2.and.phi>=xprobmin3.and.phi<=xprobmax3)weight=0.5d0
            end if
            vol=weight*ps(g)%dvolume(i,j,k);b2=sum(cart*cart)
            tot(1:4)=tot(1:4)+vol*[error**2,sum(expected*expected),error,sqrt(sum(expected*expected))]
            tot(5:7)=tot(5:7)+vol*b2*pos
            tot(8)=tot(8)+vol*b2
            if(seam)tot(9:10)=tot(9:10)+vol*[error**2,sum(expected*expected)]
            if(frame_patch==0)tot(11)=tot(11)+vol*sum(expected*expected)
            tot(12)=tot(12)+vol*sum(expected*expected)
            mins=min(mins,[ps(g)%w(i,j,k,rho_),pth(i,j,k)])
          end do
        end do
      end do
    end do
    local(9)=largest_emf
    call MPI_ALLREDUCE(local,maximum,11,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(mins,minimum,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(tot,total,12,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(stage_calls,callsmin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(stage_calls,callsmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(callsmin/=3*iit.or.callsmax/=3*iit)call mpistop('incorrect static-test RK stage count')
    if(test_world_rank==0)then
      write(history_unit,'(i6,26es25.16)')iit,qt,maximum,minimum,total
      flush(history_unit)
      if(mod(iit,50)==0.or.iit==it_max)write(*,'(a,i6,4es14.5)')'CT_TRANSPORT ',iit,qt,maximum(1),sqrt(total(1)/total(2)),minimum(2)
    end if
    if(maximum(1)>1.d-10.or.minval(minimum)<=0.d0)call mpistop('transport divergence or positivity failure')
    if(save_samples.and.any(abs(qt-[0.d0,0.404d0,0.522d0])<1.d-10))then
      write(name,'(a,i0,a,i0,a)')'sample-',iit,'-rank-',test_world_rank,'.txt'
      open(newunit=u,file=trim(name),status='replace')
      do n=1,igridstail
        g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
        do k=ixMlo3,ixMhi3
          do j=ixMlo2,ixMhi2
            do i=ixMlo1,ixMhi1
              gi=bc*bn+[i,j,k]-lo
              write(u,'(4i8,8es25.16)')frame_patch,gi,ps(g)%w(i,j,k,:)
            end do
          end do
        end do
      end do
      close(u)
    end if
    if(iit==it_max)then
      write(name,'(a,i0,a)')'transport-rank-',test_world_rank,'.txt'
      open(newunit=u,file=trim(name),status='replace')
      do n=1,igridstail
        g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
        do k=ixMlo3,ixMhi3
          do j=ixMlo2,ixMhi2
            do i=ixMlo1,ixMhi1
              gi=bc*bn+[i,j,k]-lo
              write(u,'(4i8,8es25.16)')frame_patch,gi,ps(g)%w(i,j,k,:)
            end do
          end do
        end do
      end do
      close(u)
      write(name,'(a,i0,a)')'faces-rank-',test_world_rank,'.txt'
      open(newunit=u,file=trim(name),status='replace')
      do n=1,igridstail
        g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
        do c=1,3
          fl=lo;fl(c)=fl(c)-1
          do k=fl(3),hi(3)
            do j=fl(2),hi(2)
              do i=fl(1),hi(1)
                gi=bc*bn+[i,j,k]-lo
                write(u,'(8i8,es25.16)')frame_patch,bc,c,gi,ps(g)%ws(i,j,k,c)
              end do
            end do
          end do
        end do
      end do
      close(u)
    end if
  end subroutine
end module
