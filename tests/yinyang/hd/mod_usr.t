module mod_usr
  use mod_hd
  use mod_physics, only: phys_to_conserved
  use mod_yinyang, only: yy_position, yy_basis
  implicit none
  integer :: initial_case=0
contains
  subroutine usr_init()
    call read_test_parameters()
    call set_coordinate_system('yinyang')
    call hd_activate()
    usr_init_one_grid=>initial
  end subroutine
  subroutine read_test_parameters()
    integer :: n
    namelist /yy_test/ initial_case
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,yy_test,end=10)
10    close(unitpar)
    end do
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: ix^D
    double precision :: b(3,3),v(3),bg(3),pos(3)
    v=[0.1d0,0.2d0,-0.1d0]; bg=[0.2d0,-0.1d0,0.3d0]
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call yy_basis(x(ix^D,2),x(ix^D,3),yy_patch,b)
      call yy_position(x(ix^D,1),x(ix^D,2),x(ix^D,3),yy_patch,pos)
      w(ix^D,rho_)=1.d0
      w(ix^D,p_)=1.d0
      if(initial_case==1) w(ix^D,rho_)=1.d0+0.1d0*exp(-sum((pos-[-1.2d0,0.8d0,0.d0])**2)/0.04d0)
      if(initial_case==2) then
        if(sum((pos-[-1.2d0,0.8d0,0.d0])**2)<0.04d0) w(ix^D,p_)=5.d0
      end if
      w(ix^D,mom(:))=matmul(transpose(b),v)
    {end do\}
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine
end module
