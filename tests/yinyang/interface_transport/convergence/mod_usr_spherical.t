module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved
  use mod_comm_lib, only: mpistop
  use mod_yinyang, only: yy_position, yy_basis
  implicit none
  integer :: test_case=1,direction=1,placement=1
  double precision :: pulse_width=0.18d0
  double precision :: rotation(3,3),velocity(3),origin(3)
contains
  subroutine usr_init()
    integer :: n
    double precision :: phi
    namelist /transport_test/ test_case,direction,placement,pulse_width
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,transport_test,end=10)
10    close(unitpar)
    end do
    if(abs(direction)/=1.or.test_case<0.or.test_case>1.or.placement<0.or.placement>1) &
      call mpistop('invalid transport_test configuration')
    cmax_global=0.d0 ! Before the first setdt; initial characteristic defects are zero.
    phi=0.d0
    if(placement==1) phi=146.25d0*acos(-1.d0)/180.d0
    rotation(:,1)=[-sin(phi),cos(phi),0.d0]
    rotation(:,2)=[0.d0,0.d0,1.d0]
    rotation(:,3)=[cos(phi),sin(phi),0.d0]
    origin=matmul(rotation,[-0.3d0*direction,0.d0,1.5d0])
    velocity=matmul(rotation,[dble(direction),0.d0,0.d0])
    call set_coordinate_system('spherical')
    call mhd_activate()
    usr_init_one_grid=>initial
    usr_special_bc=>boundary
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    call exact_state(0.d0,ixI^L,ixO^L,w,x)
  end subroutine
  subroutine exact_state(time,ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: time,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: ix^D
    double precision :: basis(3,3),pos(3),offset(3),field(3),phase
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call yy_basis(x(ix^D,2),x(ix^D,3),0,basis)
      call yy_position(x(ix^D,1),x(ix^D,2),x(ix^D,3),0,pos)
      offset=pos-origin-velocity*time
      w(ix^D,rho_)=1.d0
      if(test_case==0) then
        w(ix^D,rho_)=1.d0+0.5d0*exp(-sum(offset**2)/pulse_width**2)
        field=matmul(rotation,[0.2d0,0.1d0,0.3d0])
      else
        phase=0.5d0*dpi*exp(-(dot_product(offset,rotation(:,1))/pulse_width)**2)
        field=matmul(rotation,[0.d0,0.3d0*cos(phase),0.3d0*sin(phase)])
      end if
      w(ix^D,p_)=0.1d0
      w(ix^D,mom(:))=matmul(transpose(basis),velocity)
      w(ix^D,mag(:))=matmul(transpose(basis),field)
      w(ix^D,psi_)=0.d0
    {end do\}
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine
  subroutine boundary(qdt,qt,ixI^L,ixO^L,iB,w,x)
    integer,intent(in) :: ixI^L,ixO^L,iB
    double precision,intent(in) :: qdt,qt,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    call exact_state(qt,ixI^L,ixO^L,w,x)
  end subroutine
end module mod_usr
