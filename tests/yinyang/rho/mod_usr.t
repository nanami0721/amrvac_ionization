module mod_usr
  use mod_rho
  implicit none
contains
  subroutine usr_init()
    call set_coordinate_system('yinyang')
    call rho_activate()
    usr_init_one_grid=>initial
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    w(ixO^S,rho_)=1.d0
  end subroutine
end module
