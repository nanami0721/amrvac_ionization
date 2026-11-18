module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved
  use mod_yinyang, only: yy_position, yy_basis
  implicit none
contains
  subroutine usr_init()
    call set_coordinate_system('yinyang')
    call mhd_activate()
    usr_init_one_grid=>initial
    usr_special_bc=>radial_boundary
  end subroutine

  ! Normalized units: R_sun=1, polar B=1 (display scale: 10 gauss).
  ! Uniform density and pressure, no gravity/wind/heating. Analytic potential field.
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: ix^D
    double precision :: basis(3,3),position(3),normal(3),axis(3),field(3),radius
    axis=[sin(20.d0*dpi/180.d0),0.d0,cos(20.d0*dpi/180.d0)]
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      radius=x(ix^D,1)
      call yy_position(radius,x(ix^D,2),x(ix^D,3),yy_patch,position)
      call yy_basis(x(ix^D,2),x(ix^D,3),yy_patch,basis)
      normal=position/radius
      field=0.5d0*(3.d0*dot_product(axis,normal)*normal-axis)/radius**3
      w(ix^D,rho_)=1.d0
      w(ix^D,p_)=0.05d0
      w(ix^D,mom(:))=0.d0
      w(ix^D,mag(:))=matmul(transpose(basis),field)
      w(ix^D,psi_)=0.d0
    {end do\}
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine radial_boundary(qdt,qt,ixI^L,ixO^L,iB,w,x)
    integer,intent(in) :: ixI^L,ixO^L,iB
    double precision,intent(in) :: qdt,qt,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    call initial(ixI^L,ixO^L,w,x)
  end subroutine
end module
