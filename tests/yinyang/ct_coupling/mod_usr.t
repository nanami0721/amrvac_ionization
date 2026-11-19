! Exercise CT coupling workspaces on the actual Yin-Yang mesh without enabling
! unsupported MHD CT evolution. Native get_divb is used on a separate probe state.
module mod_usr
  use mod_rho
  use mod_yinyang, only: yy_ct_stage_plan,yy_prepare_ct_edges,yy_exchange_ct_edges,yy_basis
  use mod_yinyang_ct_coupling, only: yy_ct_close_ghosts,yy_ct_edge_quadrature
  use mod_functions_bfield, only: get_divb
  implicit none
  logical :: checked=.false.
contains
  subroutine usr_init()
    call set_coordinate_system('yinyang')
    call rho_activate()
    usr_init_one_grid=>initial
    usr_process_global=>check_coupling
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    w(ixO^S,rho_)=1.d0
  end subroutine
  subroutine check_coupling(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    type(yy_ct_stage_plan) :: plan
    type(ct_stage_flux),allocatable :: fluxes(:)
    type(state),target :: probe
    integer,allocatable :: receivers(:,:)
    double precision,allocatable :: saved(:,:,:,:)
    logical,allocatable :: allowed(:,:,:,:)
    integer :: n,g,i,j,k,c,axis,side,t,lo(3),hi(3),nb(3),pig(3),idx(3),q,status,jlo,jhi,pass
    double precision :: offset(3),coord(3),frame(3,3),pos(3,2),weight(3,2),exact
    double precision :: local(3),maximum(3),divb(ixG^T),dt,scale
    logical :: sides(2,2),selected,old_stagger
    if(checked) return
    checked=.true.;lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3]
    nb=[ng1(1),ng2(1),ng3(1)];dt=0.37d0;local=0.d0
    allocate(receivers(5,12*igridstail),fluxes(max_blocks));t=0
    do n=1,igridstail
      g=igrids(n);pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
      allocate(fluxes(g)%edge(ixG^T,3))
      do axis=2,3
        do side=1,2
          if(side==1.and.pig(axis)/=1) cycle
          if(side==2.and.pig(axis)/=nb(axis)) cycle
          do c=1,3
            if(c==axis) cycle ! Only tangents lie on the interface.
            idx=(lo+hi)/2
            idx(axis)=lo(axis)-1
            if(side==2) idx(axis)=hi(axis)
            t=t+1;receivers(:,t)=[g,idx,c]
          end do
        end do
      end do
    end do
    call yy_prepare_ct_edges(plan,ps,receivers(:,:t))
    if(any(.not.plan%edge%covered)) call mpistop('CT interface edge lacks donors')
    do n=1,size(plan%packed)
      g=plan%edge%samples%request(1,n);c=plan%edge%samples%request(5,n)
      idx=plan%edge%samples%request(2:4,n)+lo
      local(1)=max(local(1),abs(plan%edge%donor_length(n)-ps(g)%dsC(idx(1),idx(2),idx(3),c)))
    end do
    ! Fresh analytic stage data before each call; no interpolation feedback.
    do pass=1,3
      scale=1.d0
      if(pass==2)scale=1.7d0
      if(pass==3)scale=0.d0
      do n=1,igridstail
        g=igrids(n);fluxes(g)%edge=0.d0
        do c=1,3
          offset=0.5d0;offset(c)=0.d0
          do k=ixMlo3-1,ixMhi3
            do j=ixMlo2-1,ixMhi2
              do i=ixMlo1-1,ixMhi1
                idx=[i,j,k]
                if(idx(c)<lo(c)) cycle
                coord=ps(g)%x(i,j,k,:)+offset*ps(g)%dx(i,j,k,:)
                call yy_basis(coord(2),coord(3),yy_patch,frame)
                fluxes(g)%edge(i,j,k,c)=scale*dt*ps(g)%dsC(i,j,k,c)* &
                  dot_product([0.2d0,-0.3d0,0.5d0],frame(:,c))
              end do
            end do
          end do
        end do
      end do
      call yy_exchange_ct_edges(plan,fluxes)
      do n=1,t
        g=receivers(1,n);idx=receivers(2:4,n);c=receivers(5,n)
        offset=0.5d0;offset(c)=0.d0
        coord=ps(g)%x(idx(1),idx(2),idx(3),:)+offset*ps(g)%dx(idx(1),idx(2),idx(3),:)
        call yy_basis(coord(2),coord(3),yy_patch,frame)
        exact=scale*dt*ps(g)%dsC(idx(1),idx(2),idx(3),c)*dot_product([0.2d0,-0.3d0,0.5d0],frame(:,c))
        local(2)=max(local(2),abs(fluxes(g)%edge(idx(1),idx(2),idx(3),c)-exact))
      end do
    end do
    allocate(probe%ws(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3))
    allocate(saved(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3))
    allocate(allowed(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3))
    do n=1,igridstail
      g=igrids(n);pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
      sides(1,:)=(pig(2:3)==1);sides(2,:)=(pig(2:3)==nb(2:3))
      do c=1,3
        do k=ixGlo3-1,ixGhi3
          do j=ixGlo2-1,ixGhi2
            do i=ixGlo1-1,ixGhi1
              probe%ws(i,j,k,c)=sin(0.2d0*i+0.3d0*j-0.4d0*k+c)
            end do
          end do
        end do
      end do
      saved=probe%ws;allowed=.false.
      if(sides(1,1))allowed(ixMlo1:ixMhi1,ixGlo2-1:ixMlo2-2,ixMlo3:ixMhi3,2)=.true.
      if(sides(2,1))allowed(ixMlo1:ixMhi1,ixMhi2+1:ixGhi2,ixMlo3:ixMhi3,2)=.true.
      jlo=ixMlo2;jhi=ixMhi2
      if(sides(1,1))jlo=ixGlo2
      if(sides(2,1))jhi=ixGhi2
      if(sides(1,2))allowed(ixMlo1:ixMhi1,jlo:jhi,ixGlo3-1:ixMlo3-2,3)=.true.
      if(sides(2,2))allowed(ixMlo1:ixMhi1,jlo:jhi,ixMhi3+1:ixGhi3,3)=.true.
      do pass=1,3
        call yy_ct_close_ghosts([ixGlo1-1,ixGlo2-1,ixGlo3-1],lo,hi, &
          [ixGlo1,ixGlo2,ixGlo3],[ixGhi1,ixGhi2,ixGhi3],sides,probe%ws,ps(g)%surfaceC,status)
        if(status/=0) call mpistop('native ghost closure rejected mesh')
      end do
      if(any(probe%ws/=saved.and..not.allowed)) call mpistop('ghost closure changed protected face')
      probe%surfaceC=>ps(g)%surfaceC;probe%dvolume=>ps(g)%dvolume
      block=>probe;old_stagger=stagger_grid;stagger_grid=.true.
      call get_divb(ps(g)%w,ixG^LL,ixG^LL,divb)
      stagger_grid=old_stagger;block=>ps(g)
      do k=ixGlo3,ixGhi3
        do j=ixGlo2,ixGhi2
          selected=((j<ixMlo2.and.sides(1,1)).or.(j>ixMhi2.and.sides(2,1))) &
            .and.k>=ixMlo3.and.k<=ixMhi3
          if((k<ixMlo3.and.sides(1,2)).or.(k>ixMhi3.and.sides(2,2))) &
            selected=j>=jlo.and.j<=jhi
          if(selected)local(3)=max(local(3),maxval(abs(divb(ixMlo1:ixMhi1,j,k))))
        end do
      end do
    end do
    call MPI_ALLREDUCE(local,maximum,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(yy_world_rank==0)write(*,'(a,3es26.17)')'CT_COUPLING_MESH ',maximum
    if(maxval(maximum)>1.d-11)call mpistop('CT coupling mesh check failed')
  end subroutine
end module
