! Explicit face-fill test on a native mesh; production CT remains disabled.
module mod_usr
  use mod_rho
  use mod_yinyang, only: yy_ct_ghost_plan,yy_prepare_ct_ghosts,yy_fill_ct_ghosts,yy_basis
  use mod_functions_bfield, only: get_divb
  implicit none
  logical :: checked=.false.
contains
  subroutine usr_init()
    call set_coordinate_system('yinyang')
    call rho_activate()
    usr_init_one_grid=>initial
    usr_process_global=>check_fill
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    w(ixO^S,rho_)=1.d0
  end subroutine
  subroutine check_fill(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    type(yy_ct_ghost_plan) :: plan
    type(state),allocatable,target :: probe(:),before(:),previous(:)
    logical,allocatable :: allowed(:,:,:,:,:)
    integer :: lo(3),hi(3),fl(3),idx(3),above(3),n,g,i,j,k,c,t,p,field,closure,repeat,countall,localcount
    integer :: pig(3),nb(3),ixO^L
    double precision :: offset(3),coord(3),frame(3,3),xyz(3),vec(3),exact,error,local(3),global(3),divb(ixG^T)
    logical :: constrain,sides(2,2),old_stagger
    character(len=100) :: filename
    integer :: u
    if(checked)return
    checked=.true.;lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3];nb=[ng1(1),ng2(1),ng3(1)]
    allocate(probe(max_blocks),before(max_blocks),previous(max_blocks))
    ! Only owned blocks are populated. A sentinel poisons every unfilled halo.
    do n=1,igridstail
      g=igrids(n)
      allocate(probe(g)%ws(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3))
      probe(g)%surfaceC=>ps(g)%surfaceC;probe(g)%dvolume=>ps(g)%dvolume
    end do
    call yy_prepare_ct_ghosts(plan,ps)
    do p=1,2
      if(p==1)then
        t=size(plan%faces%radial%request,2)
      else
        t=size(plan%faces%tangent%request,2)
      end if
      do n=1,t
        if(p==1)then
          idx=plan%faces%radial%request(2:4,n)+lo;c=plan%faces%radial%request(5,n)
        else
          idx=plan%faces%tangent%request(2:4,n)+lo;c=plan%faces%tangent%request(5,n)
        end if
        fl=lo;fl(c)=fl(c)-1
        if(any(idx<fl).or.any(idx>hi))call mpistop('magnetic fill requested a ghost donor')
      end do
    end do
    allocate(allowed(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3,igridstail))
    allowed=.false.
    do n=1,size(plan%values)
      g=plan%receivers(1,n);idx=plan%receivers(2:4,n);c=plan%receivers(5,n)
      do p=1,igridstail
        if(igrids(p)==g)exit
      end do
      allowed(idx(1),idx(2),idx(3),c,p)=.true.
    end do
    write(filename,'(a,i0,a)')'ghost-rank-',yy_world_rank,'.txt'
    open(newunit=u,file=trim(filename),status='replace')
    do field=0,1
      do closure=0,1
        constrain=closure==1
        do n=1,igridstail
          g=igrids(n);probe(g)%ws=1.d20
          do c=1,3
            fl=lo;fl(c)=fl(c)-1
            do k=fl(3),hi(3)
              do j=fl(2),hi(2)
                do i=fl(1),hi(1)
                  call face_value(g,[i,j,k],c,field,exact)
                  probe(g)%ws(i,j,k,c)=exact
                end do
              end do
            end do
          end do
          before(g)%ws=probe(g)%ws
        end do
        do repeat=1,3
          call yy_fill_ct_ghosts(plan,probe,constrain)
          do n=1,igridstail
            g=igrids(n)
            if(any(probe(g)%ws/=before(g)%ws.and..not.allowed(:,:,:,:,n))) &
              call mpistop('magnetic fill changed a protected face')
            if(repeat==1)then
              previous(g)%ws=probe(g)%ws
            else
              if(any(probe(g)%ws/=previous(g)%ws))call mpistop('magnetic fill not reproducible')
            end if
          end do
        end do
        ! Keep the fully allocated test state alive across several watchdog
        ! samples; otherwise sub-second MPI jobs can hide their peak residency.
        if(field==0.and.closure==0) call sleep(1)
        local=0.d0;localcount=size(plan%values)
        do n=1,size(plan%values)
          g=plan%receivers(1,n);idx=plan%receivers(2:4,n);c=plan%receivers(5,n)
          call face_value(g,idx,c,field,exact)
          error=abs(probe(g)%ws(idx(1),idx(2),idx(3),c)-exact)
          local(1)=max(local(1),error);local(2)=local(2)+error**2
          ! Global logical face IDs make the output independent of MPI ownership.
          idx=idx-lo+[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]* &
            [block_nx1,block_nx2,block_nx3]
          write(u,'(7i8,es26.17)')field,closure,yy_patch,idx,c,probe(g)%ws( &
            plan%receivers(2,n),plan%receivers(3,n),plan%receivers(4,n),c)
        end do
        do n=1,igridstail
          g=igrids(n);pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
          sides(1,:)=(pig(2:3)==1);sides(2,:)=(pig(2:3)==nb(2:3))
          block=>probe(g);old_stagger=stagger_grid;stagger_grid=.true.
          ! Only evaluate populated slabs, not unrelated sentinel halos.
          do p=1,2
            if(sides(p,1))then
              ixO^L=ixM^LL;
              ixOmin2=ixGlo2;ixOmax2=ixMlo2-1
              if(p==2)then
                ixOmin2=ixMhi2+1;ixOmax2=ixGhi2
              end if
              call get_divb(ps(g)%w,ixG^LL,ixO^L,divb)
              local(3)=max(local(3),maxval(abs(divb(ixO^S))))
            end if
            if(sides(p,2))then
              ixO^L=ixM^LL;
              if(sides(1,1))ixOmin2=ixGlo2
              if(sides(2,1))ixOmax2=ixGhi2
              ixOmin3=ixGlo3;ixOmax3=ixMlo3-1
              if(p==2)then
                ixOmin3=ixMhi3+1;ixOmax3=ixGhi3
              end if
              call get_divb(ps(g)%w,ixG^LL,ixO^L,divb)
              local(3)=max(local(3),maxval(abs(divb(ixO^S))))
            end if
          end do
          stagger_grid=old_stagger;block=>ps(g)
        end do
        call MPI_ALLREDUCE(local(1),global(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
        call MPI_ALLREDUCE(local(2),global(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
        call MPI_ALLREDUCE(local(3),global(3),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
        call MPI_ALLREDUCE(localcount,countall,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
        global(2)=sqrt(global(2)/countall)
        if(yy_world_rank==0)write(*,'(a,3i9,3es26.17)')'MAGNETIC_GHOST ',field,closure,countall,global
        if(global(1)>0.2d0)call mpistop('magnetic ghost field inaccurate or used poisoned donors')
        if(constrain.and.global(3)>1.d-11)call mpistop('magnetic ghost divergence constraint failed')
      end do
    end do
    close(u)
  contains
    subroutine face_value(g,idx,c,field,value)
      integer,intent(in)::g,idx(3),c,field
      double precision,intent(out)::value
      integer::above(3)
      double precision::offset(3),coord(3),frame(3,3),xyz(3),vec(3)
      above=max(idx,[ixGlo1,ixGlo2,ixGlo3]);offset=0.d0;offset(c)=0.5d0
      if(idx(c)<above(c))offset(c)=-0.5d0
      coord=ps(g)%x(above(1),above(2),above(3),:)+offset*ps(g)%dx(above(1),above(2),above(3),:)
      call yy_basis(coord(2),coord(3),yy_patch,frame);xyz=coord(1)*frame(:,1)
      vec=[0.2d0,-0.3d0,0.5d0]
      if(field==1)vec=[xyz(2),xyz(3),xyz(1)]
      value=dot_product(vec,frame(:,c))
    end subroutine
  end subroutine
end module
