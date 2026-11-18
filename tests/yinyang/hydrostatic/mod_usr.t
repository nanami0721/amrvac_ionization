module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved,phys_to_primitive
  use mod_comm_lib, only: mpistop
  use mod_yinyang, only: yy_position,yy_basis
  implicit none
  integer :: equilibrium_case=2
  logical :: spherical_control=.false.,diagnostics=.true.,diagnostic_every_step=.false.
  double precision,parameter :: gm=1.d0,ct2=0.5d0,rho_base=1.d0
  double precision :: axis(3),tau,last_dump=-huge(1.d0)
  double precision :: beta_reference=4.d0,bpole,output_interval=1.d0
  double precision,allocatable :: weighted_volume(:,:,:,:),budget(:,:),before_source(:,:)
  integer :: budget_unit,rank_id
contains
  subroutine usr_init()
    integer :: n
    namelist /equilibrium_test/ equilibrium_case,spherical_control,diagnostics,diagnostic_every_step, &
      beta_reference,output_interval
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,equilibrium_test,end=10)
10    close(unitpar)
    end do
    if(equilibrium_case<0.or.equilibrium_case>3) call mpistop('invalid equilibrium_case')
    if(.not.(beta_reference>0.d0).or..not.(output_interval>0.d0)) &
      call mpistop('beta and diagnostic interval must be positive')
    bpole=sqrt(2.d0*ct2*rho_base/beta_reference)
    axis=[sin(20.d0*acos(-1.d0)/180.d0),0.d0,cos(20.d0*acos(-1.d0)/180.d0)]
    if(spherical_control) then
      call set_coordinate_system('spherical')
    else
      call set_coordinate_system('yinyang')
    end if
    usr_gravity=>central_gravity
    call mhd_activate()
    usr_set_parameters=>parameters
    usr_init_one_grid=>initial
    usr_special_bc=>boundary
    if(diagnostics) then
      usr_process_global=>begin_step
      usr_process_adv_global=>end_step
      usr_fv_stage_diagnostic=>stage_budget
    end if
  end subroutine
  subroutine parameters()
    if(any(boundary_divbfix)) call mpistop('equilibrium test prescribes analytic B; disable boundary_divbfix')
    if(mhd_gravity.neqv.(equilibrium_case==0.or.equilibrium_case==2)) call mpistop('gravity flag disagrees with equilibrium_case')
    if(t_stepper/=threestep.or.t_integrator/=ssprk3.or.ssprk_order/=3.or. &
       dimsplit.or.local_timestep) call mpistop('equilibrium budget requires unsplit SSPRK(3,3)')
    tau=1.d0/sqrt(eos%gamma*ct2)
  end subroutine
  subroutine central_gravity(ixI^L,ixO^L,w,x,gravity_field)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: w(ixI^S,nw),x(ixI^S,ndim)
    double precision,intent(out) :: gravity_field(ixI^S,ndim)
    gravity_field(ixO^S,1)=-gm/x(ixO^S,1)**2
    gravity_field(ixO^S,2)=0.d0
    gravity_field(ixO^S,3)=0.d0
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: ix^D,patch
    double precision :: b(3,3),pos(3),normal(3),field(3),rad
    patch=0
    if(yinyang_enabled) patch=yy_patch
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      rad=x(ix^D,1)
      call yy_position(rad,x(ix^D,2),x(ix^D,3),patch,pos)
      call yy_basis(x(ix^D,2),x(ix^D,3),patch,b)
      normal=pos/rad
      w(ix^D,rho_)=rho_base
      if(equilibrium_case==0.or.equilibrium_case==2) w(ix^D,rho_)=rho_base*exp(gm/ct2*(1.d0/rad-1.d0))
      w(ix^D,p_)=ct2*w(ix^D,rho_)
      w(ix^D,mom(:))=0.d0
      field=0.d0
      if(equilibrium_case/=0) field=0.5d0*bpole*(3.d0*dot_product(axis,normal)*normal-axis)/rad**3
      if(equilibrium_case==3) field=bpole*axis
      w(ix^D,mag(:))=matmul(transpose(b),field)
      w(ix^D,psi_)=0.d0
    {end do\}
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine
  subroutine boundary(qdt,qt,ixI^L,ixO^L,iB,w,x)
    integer,intent(in) :: ixI^L,ixO^L,iB
    double precision,intent(in) :: qdt,qt,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    call initial(ixI^L,ixO^L,w,x)
  end subroutine
  subroutine initialize_diagnostics()
    integer :: n,g,i,j,k,a,c
    double precision :: dteta,dphi,tt,pp,ot,op,solid,weight
    character(len=80) :: filename
    allocate(weighted_volume(ixG^T,maxval(igrids(1:igridstail))))
    allocate(budget(6,size(weighted_volume,4)),before_source(2,size(weighted_volume,4)))
    budget=0.d0;before_source=0.d0
    do n=1,igridstail
      g=igrids(n)
      do k=ixMlo3,ixMhi3
        do j=ixMlo2,ixMhi2
          weight=1.d0
          if(yinyang_enabled) then
            dteta=rnode(rpdx2_,g);dphi=rnode(rpdx3_,g);weight=0.d0
            solid=2.d0*sin(ps(g)%x(ixMlo1,j,k,2))*sin(dteta/2)*dphi
            do a=0,15
              tt=ps(g)%x(ixMlo1,j,k,2)+((a+0.5d0)/16-0.5d0)*dteta
              do c=0,15
                pp=ps(g)%x(ixMlo1,j,k,3)+((c+0.5d0)/16-0.5d0)*dphi
                ot=atan2(sqrt((sin(tt)*cos(pp))**2+cos(tt)**2),sin(tt)*sin(pp))
                op=atan2(cos(tt),-sin(tt)*cos(pp))
                if(ot>=xprobmin2.and.ot<=xprobmax2.and.op>=xprobmin3.and.op<=xprobmax3) then
                  weight=weight+0.5d0*(cos(tt-dteta/32)-cos(tt+dteta/32))*dphi/16
                else
                  weight=weight+(cos(tt-dteta/32)-cos(tt+dteta/32))*dphi/16
                end if
              end do
            end do
            weight=weight/solid
          end if
          do i=ixMlo1,ixMhi1
            weighted_volume(i,j,k,g)=weight*ps(g)%dvolume(i,j,k)
          end do
        end do
      end do
    end do
    rank_id=mype
    if(yinyang_enabled) rank_id=yy_world_rank
    write(filename,'(a,i4.4,a)') 'budget_r',rank_id,'.txt'
    open(newunit=budget_unit,file=trim(filename),status='new')
    write(budget_unit,'(a)') '# time mass energy radial_mass radial_energy gravity_work '// &
      'all_source_energy weighted_flux_mass weighted_flux_energy'
    call verify_gravity_source()
  end subroutine
  subroutine verify_gravity_source()
    use mod_gravity, only: gravity_add_source
    integer :: g,ix^D,patch
    double precision :: w0(ixG^T,nw),wp(ixG^T,nw),w1(ixG^T,nw),expected,err
    double precision :: gravity_field(ixG^T,3),basis(3,3),position(3),cart_gravity(3)
    logical :: active
    if(.not.mhd_gravity) return
    g=igrids(1)
    call initial(ixG^LL,ixG^LL,w0,ps(g)%x)
    wp=w0
    call phys_to_primitive(ixG^LL,ixG^LL,wp,ps(g)%x)
    wp(ixG^T,mom(1))=0.2d0;wp(ixG^T,mom(2))=-0.3d0;wp(ixG^T,mom(3))=0.4d0
    w0=wp
    call phys_to_conserved(ixG^LL,ixG^LL,w0,ps(g)%x)
    w1=w0;active=.false.;err=0.d0
    patch=0
    if(yinyang_enabled) patch=yy_patch
    call central_gravity(ixG^LL,ixM^LL,w0,ps(g)%x,gravity_field)
    call gravity_add_source(0.01d0,ixG^LL,ixM^LL,w0,wp,w1,ps(g)%x,.true.,.false.,active)
    {do ix^DB=ixMlo^DB,ixMhi^DB\}
      call yy_basis(ps(g)%x(ix^D,2),ps(g)%x(ix^D,3),patch,basis)
      call yy_position(ps(g)%x(ix^D,1),ps(g)%x(ix^D,2),ps(g)%x(ix^D,3),patch,position)
      cart_gravity=matmul(basis,gravity_field(ix^D,:))
      err=max(err,maxval(abs(cart_gravity+gm*position/ps(g)%x(ix^D,1)**3)))
      err=max(err,abs(w1(ix^D,mom(2))-w0(ix^D,mom(2))),abs(w1(ix^D,mom(3))-w0(ix^D,mom(3))))
      expected=-0.01d0*gm*w0(ix^D,rho_)/ps(g)%x(ix^D,1)**2
      err=max(err,abs(w1(ix^D,mom(1))-w0(ix^D,mom(1))-expected))
      err=max(err,abs(w1(ix^D,e_)-w0(ix^D,e_)-0.2d0*expected))
    {end do\}
    if(.not.active.or.err>1.d-12) call mpistop('central gravity source self-check failed')
    if(rank_id==0) write(*,'(a,es12.4)') 'central gravity coordinate/source self-check error: ',err
  end subroutine
  subroutine stage_budget(phase,qdt,ixI^L,ixO^L,idims^LIM,fC,wCT,w,x)
    integer,intent(in) :: phase,ixI^L,ixO^L,idims^LIM
    double precision,intent(in) :: qdt,fC(ixI^S,nwflux,ndim),wCT(ixI^S,nw),w(ixI^S,nw),x(ixI^S,ndim)
    integer :: g,ix^D,iw,component
    double precision :: factor,vol,part,total(2),radial(2),work,fluxdiv
    if(.not.allocated(budget)) call mpistop('budget must be initialized before finite-volume stages')
    if(idimsmin/=1.or.idimsmax/=3) call mpistop('budget needs all spatial directions')
    select case(istep)
    case(1)
      factor=1.d0/6.d0
    case(2)
      factor=2.d0/3.d0
    case(3)
      factor=1.d0
    case default
      call mpistop('unsupported stage in equilibrium budget')
    end select
    g=block%igrid
    total(1)=sum(w(ixO^S,rho_)*weighted_volume(ixO^S,g))
    total(2)=sum(w(ixO^S,e_)*weighted_volume(ixO^S,g))
    if(phase==fv_after_sources) then
      budget(4,g)=budget(4,g)+factor*(total(2)-before_source(2,g))
      return
    end if
    before_source(:,g)=total;radial=0.d0;work=0.d0;total=0.d0
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      vol=weighted_volume(ix^D,g);part=vol/block%dvolume(ix^D)
      if(mhd_gravity) work=work-qdt*gm*wCT(ix^D,mom(1))*vol/x(ix^D,1)**2
      do component=1,2
        iw=rho_
        if(component==2) iw=e_
        fluxdiv=fC(ix1,ix2,ix3,iw,1)-fC(ix1-1,ix2,ix3,iw,1) &
               +fC(ix1,ix2,ix3,iw,2)-fC(ix1,ix2-1,ix3,iw,2) &
               +fC(ix1,ix2,ix3,iw,3)-fC(ix1,ix2,ix3-1,iw,3)
        total(component)=total(component)+part*fluxdiv
        if(abs(x(ix^D,1)-0.5d0*dxlevel(1)-xprobmin1)<1.d-10) &
          radial(component)=radial(component)-part*fC(ix1-1,ix2,ix3,iw,1)
        if(abs(x(ix^D,1)+0.5d0*dxlevel(1)-xprobmax1)<1.d-10) &
          radial(component)=radial(component)+part*fC(ix1,ix2,ix3,iw,1)
      end do
    {end do\}
    budget(1:2,g)=budget(1:2,g)+factor*radial
    budget(3,g)=budget(3,g)+factor*work
    budget(5:6,g)=budget(5:6,g)+factor*total
  end subroutine
  subroutine begin_step(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    if(.not.allocated(budget)) call initialize_diagnostics()
    call record_state(iit,qt)
  end subroutine
  subroutine end_step(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    call record_state(iit+1,qt+dt)
  end subroutine
  subroutine record_state(iit,qt)
    use mod_geometry, only: curlvector
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    integer :: n,g,u,k,ix^D,header(5),idirmin,patch
    double precision :: mass,energy,sums(6),old_dx(ndim),divb(ixG^T),cur(ixG^T,3),bfield(ixG^T,3)
    double precision :: reference(ixG^T,nw),divref(ixG^T),curref(ixG^T,3)
    double precision :: values(12,block_nx1*block_nx2*block_nx3)
    type(state),pointer :: old_block
    character(len=100) :: filename
    if(abs(qt-last_dump)<1.d-12) return
    if(.not.diagnostic_every_step.and.abs(qt/(tau*output_interval)-nint(qt/(tau*output_interval)))>1.d-10) return
    last_dump=qt;mass=0.d0;energy=0.d0;sums=0.d0
    do n=1,igridstail
      g=igrids(n)
      mass=mass+sum(ps(g)%w(ixM^T,rho_)*weighted_volume(ixM^T,g))
      energy=energy+sum(ps(g)%w(ixM^T,e_)*weighted_volume(ixM^T,g))
      sums=sums+budget(:,g)
    end do
    write(budget_unit,'(9es25.16)') qt,mass,energy,sums
    flush(budget_unit)
    patch=0
    if(yinyang_enabled) patch=yy_patch
    write(filename,'(a,i1,a,i4.4,a,i7.7,a)') 'equilibrium_p',patch,'_r',rank_id,'_i',iit,'.bin'
    open(newunit=u,file=trim(filename),access='stream',form='unformatted',status='new')
    header=[1,patch,rank_id,iit,igridstail*block_nx1*block_nx2*block_nx3]
    write(u) header,qt
    old_block=>block;old_dx=dxlevel
    do n=1,igridstail
      g=igrids(n);block=>ps(g)
      ^D&dxlevel(^D)=rnode(rpdx^D_,g);
      call get_divb(ps(g)%w,ixG^LL,ixM^LL,divb)
      bfield=ps(g)%w(ixG^T,mag(:))
      call curlvector(bfield,ixG^LL,ixM^LL,cur,idirmin,1,3)
      call initial(ixG^LL,ixG^LL,reference,ps(g)%x)
      call get_divb(reference,ixG^LL,ixM^LL,divref)
      bfield=reference(ixG^T,mag(:))
      call curlvector(bfield,ixG^LL,ixM^LL,curref,idirmin,1,3)
      k=0
      {do ix^DB=ixMlo^DB,ixMhi^DB\}
        k=k+1
        values(:,k)=[ps(g)%x(ix^D,1),ps(g)%x(ix^D,2),ps(g)%x(ix^D,3), &
          divb(ix^D),cur(ix^D,1),cur(ix^D,2),cur(ix^D,3), &
          divref(ix^D),curref(ix^D,1),curref(ix^D,2),curref(ix^D,3),weighted_volume(ix^D,g)]
      {end do\}
      write(u) values
    end do
    close(u);block=>old_block;dxlevel=old_dx
  end subroutine
end module mod_usr
