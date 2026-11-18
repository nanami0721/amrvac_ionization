! Native spherical CT prerequisite; this does not enable Yin-Yang CT.
module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved
  use mod_yinyang, only: yy_basis,yy_position
  use mod_functions_bfield, only: get_divb
  use mod_comm_lib, only: mpistop
  implicit none
  integer :: stage_calls=0,ct_test_case=0
  double precision :: emf_scale=1.d0
  logical :: scale_in_stage=.false.
contains
  subroutine usr_init()
    integer :: n
    namelist /ct_test/ ct_test_case,emf_scale,scale_in_stage
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,ct_test,end=10)
10    close(unitpar)
    end do
    if(ct_test_case<0.or.ct_test_case>1) call mpistop('invalid CT test case')
    usr_set_parameters=>validate_stage
    if(.not.scale_in_stage) usr_set_electric_field=>scale_edges
    call set_coordinate_system('spherical')
    call mhd_activate()
    if(.not.stagger_grid) call mpistop('CT preflight requires staggered fields')
    usr_init_one_grid=>initial
    usr_init_vector_potential=>potential
    usr_ct_stage=>check_stage
    usr_process_global=>diagnostic
    usr_process_adv_global=>diagnostic
  end subroutine
  subroutine validate_stage()
    if(scale_in_stage.and..not.ct_stage_split) call mpistop('stage scaling requires split CT')
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: ix^D
    double precision :: frame(3,3)
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    w(ixO^S,mom(:))=0.d0
    if(ct_test_case==1) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        call yy_basis(x(ix^D,2),x(ix^D,3),0,frame)
        w(ix^D,mom(:))=matmul(transpose(frame),[0.05d0,0.02d0,-0.03d0])
      {end do\}
    end if
    call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
    call mhd_face_to_center(ixO^L,block)
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine
  subroutine potential(ixI^L,ixO^L,x,A,idir)
    integer,intent(in) :: ixI^L,ixO^L,idir
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(out) :: A(ixI^S)
    integer :: ix^D
    double precision :: frame(3,3),pos(3),cart(3)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call yy_basis(x(ix^D,2),x(ix^D,3),0,frame)
      call yy_position(x(ix^D,1),x(ix^D,2),x(ix^D,3),0,pos)
      ! A = B cross r / 2, B=(0.2,-0.1,0.3).
      cart=0.5d0*[-0.1d0*pos(3)-0.3d0*pos(2), &
        0.3d0*pos(1)-0.2d0*pos(3),0.2d0*pos(2)+0.1d0*pos(1)]
      if(ct_test_case==1) cart(3)=cart(3)+0.02d0*exp(-sum((pos-[1.4d0,0.d0,0.d0])**2)/0.2d0**2)
      A(ix^D)=dot_product(cart,frame(:,idir))
    {end do\}
  end subroutine
  subroutine scale_edges(ixI^L,ixO^L,qt,qdt,fE,s)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: qt,qdt
    double precision,intent(inout) :: fE(ixI^S,sdim:3)
    type(state) :: s
    integer :: idir,ixC^L
    do idir=sdim,3
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);
      fE(ixC^S,idir)=emf_scale*fE(ixC^S,idir)
    end do
  end subroutine
  subroutine check_stage(qt,qdt,sCT,s,fluxes)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    double precision,intent(in) :: qt,qdt
    type(state),intent(in) :: sCT(max_blocks),s(max_blocks)
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    integer :: n,g,lo,hi
    stage_calls=stage_calls+1
    call MPI_ALLREDUCE(stage_calls,lo,1,MPI_INTEGER,MPI_MIN,icomm,ierrmpi)
    call MPI_ALLREDUCE(stage_calls,hi,1,MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
    if(lo/=hi) call mpistop('CT stage callback rank counts disagree')
    do n=1,igridstail_active
      g=igrids_active(n)
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT stage edge workspace missing')
      if(.not.all(ieee_is_finite(fluxes(g)%edge))) call mpistop('CT stage edge workspace uninitialized')
      if(.not.all(ieee_is_finite(fluxes(g)%face))) call mpistop('CT stage face workspace uninitialized')
      if(scale_in_stage) fluxes(g)%edge=emf_scale*fluxes(g)%edge
    end do
    if(mype==0) write(*,'(a,i8)') 'CT_STAGE_READY ',stage_calls
  end subroutine
  subroutine diagnostic(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    double precision :: divb(ixG^T),local,maximum
    integer :: n,g
    local=0.d0
    do n=1,igridstail
      g=igrids(n);block=>ps(g)
      call get_divb(ps(g)%w,ixG^LL,ixM^LL,divb)
      local=max(local,maxval(abs(divb(ixM^T))))
    end do
    call MPI_ALLREDUCE(local,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
    if(mype==0) write(*,'(a,i8,2es24.15)') 'CT_PREFLIGHT ',iit,qt,maximum
    if(maximum>1.d-10) call mpistop('native spherical CT divergence invariant failed')
  end subroutine
end module
