module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved
  use mod_comm_lib, only: mpistop
  use mod_yinyang, only: yy_position, yy_basis
  implicit none
  integer :: test_case=0,direction=1,placement=1,radial_bc=0
  logical :: native_diagnostics=.true.,native_every_step=.false.
  double precision,allocatable :: native_saved(:,:,:,:,:)
  integer :: last_native_it=-1
  double precision :: rotation(3,3),velocity(3),origin(3)
contains
  subroutine usr_init()
    integer :: n
    double precision :: phi
    namelist /transport_test/ test_case,direction,placement,radial_bc, &
      native_diagnostics,native_every_step
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,transport_test,end=10)
10    close(unitpar)
    end do
    if(abs(direction)/=1.or.test_case<0.or.test_case>1.or.placement<0.or.placement>1) &
      call mpistop('invalid transport_test configuration')
    if(radial_bc<0.or.radial_bc>1) call mpistop('invalid radial_bc')
    cmax_global=0.d0 ! Before the first setdt; initial characteristic defects are zero.
    phi=0.d0
    if(placement==1) phi=146.25d0*acos(-1.d0)/180.d0
    rotation(:,1)=[-sin(phi),cos(phi),0.d0]
    rotation(:,2)=[0.d0,0.d0,1.d0]
    rotation(:,3)=[cos(phi),sin(phi),0.d0]
    origin=matmul(rotation,[-0.3d0*direction,0.d0,1.5d0])
    velocity=matmul(rotation,[dble(direction),0.d0,0.d0])
    call set_coordinate_system('yinyang')
    call mhd_activate()
    usr_init_one_grid=>initial
    usr_special_bc=>boundary
    if(native_diagnostics) then
      usr_process_global=>native_begin
      usr_process_adv_global=>native_end
    end if
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
      call yy_basis(x(ix^D,2),x(ix^D,3),yy_patch,basis)
      call yy_position(x(ix^D,1),x(ix^D,2),x(ix^D,3),yy_patch,pos)
      offset=pos-origin-velocity*time
      w(ix^D,rho_)=1.d0
      if(test_case==0) then
        w(ix^D,rho_)=1.d0+0.5d0*exp(-sum(offset**2)/0.18d0**2)
        field=matmul(rotation,[0.2d0,0.1d0,0.3d0])
      else
        phase=0.5d0*dpi*exp(-(dot_product(offset,rotation(:,1))/0.18d0)**2)
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
    integer :: ix^D,ir
    double precision :: sign_out
    double precision :: b_reference,b_old,b_new,outgoing
    if(radial_bc/=0.and.any(boundary_divbfix(1:2))) &
      call mpistop('characteristic Bn and radial divB closure cannot be combined')
    call exact_state(qt,ixI^L,ixO^L,w,x)
    if(radial_bc==0.or.cmax_global<=0.d0) return
    select case(iB)
    case(1)
      ir=ixOmax1+1;sign_out=-1.d0
    case(2)
      ir=ixOmin1-1;sign_out=1.d0
    case default
      call mpistop('radial GLM boundary called on angular face')
    end select
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! These angular corner ghosts are overwritten by overset exchange.
      ! Do not read their possibly stale interior-radial reference values.
      if(x(ix^D,2)<xprobmin2.or.x(ix^D,2)>xprobmax2.or. &
         x(ix^D,3)<xprobmin3.or.x(ix^D,3)>xprobmax3) cycle
      b_reference=reference_br(qt,x(ir,ix2,ix3,1),x(ir,ix2,ix3,2),x(ir,ix2,ix3,3))
      ! Outgoing defect in outward-normal coordinates, incoming defect = 0.
      outgoing=w(ir,ix2,ix3,psi_)+cmax_global*sign_out*(w(ir,ix2,ix3,mag(1))-b_reference)
      b_old=w(ix^D,mag(1))
      b_new=b_old+sign_out*outgoing/(2.d0*cmax_global)
      w(ix^D,mag(1))=b_new
      w(ix^D,psi_)=0.5d0*outgoing
      ! Preserve the imposed thermal pressure when changing magnetic energy.
      w(ix^D,e_)=w(ix^D,e_)+0.5d0*(b_new**2-b_old**2)
    {end do\}
  end subroutine

  pure function reference_br(qt,r,theta,phi) result(br)
    double precision,intent(in) :: qt,r,theta,phi
    double precision :: br,basis(3,3),pos(3),offset(3),field(3),phase
    call yy_basis(theta,phi,yy_patch,basis)
    call yy_position(r,theta,phi,yy_patch,pos)
    offset=pos-origin-velocity*qt
    if(test_case==0) then
      field=matmul(rotation,[0.2d0,0.1d0,0.3d0])
    else
      phase=0.5d0*dpi*exp(-(dot_product(offset,rotation(:,1))/0.18d0)**2)
      field=matmul(rotation,[0.d0,0.3d0*cos(phase),0.3d0*sin(phase)])
    end if
    br=dot_product(basis(:,1),field)
  end function

  subroutine native_begin(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    if(iit==0) call native_dump(iit,qt)
  end subroutine
  subroutine native_end(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    call native_dump(iit+1,qt+dt)
  end subroutine
  subroutine native_dump(iit,qt)
    use mod_ghostcells_update, only: getbc
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    integer :: g,n,u,ios,ix^D,k,header(5)
    double precision :: dnum(ixG^T),dexact(ixG^T),reference(ixG^T,nw)
    double precision :: values(10,block_nx1*block_nx2*block_nx3)
    double precision :: old_dx(ndim)
    type(state),pointer :: old_block
    character(len=128) :: filename
    if(iit==last_native_it) return
    if(.not.native_every_step.and.abs(qt-0.3d0*nint(qt/0.3d0))>1.d-10) return
    last_native_it=iit
    old_block=>block;old_dx=dxlevel
    ! Diagnostic ghost refresh must not change the simulation, especially when
    ! radial divergence closure depends on the previous angular corner ghosts.
    if(.not.allocated(native_saved)) allocate(native_saved(ixG^T,nw,igridstail))
    if(size(native_saved,5)/=igridstail) call mpistop('native diagnostic requires fixed block ownership')
    do n=1,igridstail
      native_saved(:,:,:,:,n)=ps(igrids(n))%w
    end do
    call getbc(qt,0.d0,ps,iwstart,nwgc)
    write(filename,'(a,i1,a,i4.4,a,i7.7,a)') 'native_p',yy_patch,'_r',yy_world_rank,'_i',iit,'.bin'
    open(newunit=u,file=trim(filename),status='new',access='stream',form='unformatted',iostat=ios)
    if(ios/=0) call mpistop('cannot create native diagnostic file')
    header=[1,yy_patch,yy_world_rank,iit,igridstail*block_nx1*block_nx2*block_nx3]
    write(u) header,qt,cmax_global
    do n=1,igridstail
      g=igrids(n);block=>ps(g)
      ^D&dxlevel(^D)=rnode(rpdx^D_,g);
      call get_divb(ps(g)%w,ixG^LL,ixM^LL,dnum)
      call exact_state(qt,ixG^LL,ixG^LL,reference,ps(g)%x)
      call get_divb(reference,ixG^LL,ixM^LL,dexact)
      k=0
      {do ix^DB=ixMlo^DB,ixMhi^DB\}
        k=k+1
        values(:,k)=[ps(g)%x(ix^D,1),ps(g)%x(ix^D,2),ps(g)%x(ix^D,3), &
          dnum(ix^D),dexact(ix^D),block%dvolume(ix^D),ps(g)%w(ix^D,psi_), &
          ps(g)%w(ix^D,mag(1))-reference(ix^D,mag(1)), &
          ps(g)%w(ix^D,mag(2))-reference(ix^D,mag(2)), &
          ps(g)%w(ix^D,mag(3))-reference(ix^D,mag(3))]
      {end do\}
      write(u) values
    end do
    close(u)
    do n=1,igridstail
      ps(igrids(n))%w=native_saved(:,:,:,:,n)
    end do
    block=>old_block;dxlevel=old_dx
  end subroutine
end module
