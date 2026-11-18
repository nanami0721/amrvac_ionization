module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved
  use mod_comm_lib, only: mpistop
  implicit none
  integer :: direction=1
  double precision :: pulse_width=0.18d0
contains
  subroutine usr_init()
    integer :: n
    namelist /transport_test/ direction,pulse_width
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,transport_test,end=10)
10    close(unitpar)
    end do
    if(abs(direction)/=1.or.pulse_width<=0.d0) call mpistop('invalid transport parameters')
    call set_coordinate_system('Cartesian_1.75D')
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
    double precision :: phase(ixI^S)
    phase(ixO^S)=0.5d0*dpi*exp(-((x(ixO^S,1)+0.3d0*direction-direction*time)/pulse_width)**2)
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=0.1d0
    w(ixO^S,mom(1))=dble(direction)
    w(ixO^S,mom(2))=0.d0
    w(ixO^S,mom(3))=0.d0
    w(ixO^S,mag(1))=0.d0
    w(ixO^S,mag(2))=0.3d0*cos(phase(ixO^S))
    w(ixO^S,mag(3))=0.3d0*sin(phase(ixO^S))
    if(psi_>0) w(ixO^S,psi_)=0.d0
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine
  subroutine boundary(qdt,qt,ixI^L,ixO^L,iB,w,x)
    integer,intent(in) :: ixI^L,ixO^L,iB
    double precision,intent(in) :: qdt,qt,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    call exact_state(qt,ixI^L,ixO^L,w,x)
  end subroutine
end module mod_usr
