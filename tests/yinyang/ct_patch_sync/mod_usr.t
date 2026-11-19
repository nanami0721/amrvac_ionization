module mod_usr
  use mod_rho
  use mod_yinyang, only: yy_ct_patch_plan,yy_prepare_ct_patch,yy_sync_ct_edges,yy_ordered_ct_edges, &
    yy_ct_stage_plan,yy_prepare_ct_edges,yy_ct_edge_is_owner,yy_ct_ghost_plan,yy_prepare_ct_ghosts, &
    yy_fill_ct_ghosts,yy_sync_ct_faces,yy_ordered_ct_ghosts,yy_basis
  use mod_yinyang_patch_sync, only: yy_patch_copy_exchange
  use mod_yinyang_ct_coupling, only: yy_ct_edge_quadrature
  use mod_functions_bfield, only: get_divb
  implicit none
  logical  ::  checked=.false.
contains
  subroutine usr_init()
    call set_coordinate_system('yinyang');call rho_activate()
    usr_init_one_grid=>initial;usr_process_global=>check_sync
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    w(ixO^S,rho_)=1.d0
  end subroutine
  subroutine check_sync(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    type(yy_ct_patch_plan) :: patch
    type(yy_ct_stage_plan) :: cross
    type(yy_ct_ghost_plan) :: ghost
    type(ct_stage_flux),allocatable :: fluxes(:),saved(:)
    type(state),allocatable,target :: probe(:),base(:),manual(:)
    integer,allocatable :: receivers(:,:)
    integer :: lo(3),hi(3),bn(3),nn(3),nb(3),bc(3),gi(3),idx(3),fl(3),owner(3)
    integer :: n,g,i,j,k,c,d,t,pass,field,u,bad,allbad,holes,allholes,mixed,allmixed,corners,allcorners
    double precision :: offset(3),pos(3),frame(3,3),quad(3,2),weights(3,2),value,local(2),maximum(2),divb(ixG^T)
    logical :: is_interface,old_stagger,physical,active
    character(len=100) :: filename
    if(checked)return
    checked=.true.;lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3]
    bn=[block_nx1,block_nx2,block_nx3];nn=[domain_nx1,domain_nx2,domain_nx3];nb=nn/bn
    allocate(fluxes(max_blocks),saved(max_blocks),probe(max_blocks),base(max_blocks),manual(max_blocks))
    call yy_prepare_ct_patch(patch,ps)
    call yy_prepare_ct_ghosts(ghost,ps,face_moments=.true.)
    bad=0;local=0.d0;corners=0
    do n=1,igridstail
      g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
      allocate(fluxes(g)%edge(ixG^T,3),fluxes(g)%face(ixG^T,1,3))
      fluxes(g)%edge=-1.d20;fluxes(g)%face=123.d0
      do c=1,3
        fl=lo-1;fl(c)=lo(c)
        do k=fl(3),hi(3)
          do j=fl(2),hi(2)
            do i=fl(1),hi(1)
              gi=bc*bn+[i,j,k]-lo
              fluxes(g)%edge(i,j,k,c)=coded_edge(gi,c,bc)
            end do
          end do
        end do
      end do
    end do
    call yy_sync_ct_edges(patch,fluxes)
    do n=1,igridstail
      g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
      do c=1,3
        fl=lo-1;fl(c)=lo(c)
        do k=ixGlo3,ixGhi3
          do j=ixGlo2,ixGhi2
            do i=ixGlo1,ixGhi1
              idx=[i,j,k]
              if(any(idx<fl).or.any(idx>hi))then
                if(fluxes(g)%edge(i,j,k,c)/=-1.d20)bad=bad+1
                cycle
              end if
              gi=bc*bn+idx-lo;owner=max(0,min(nb-1,gi/bn))
              if(fluxes(g)%edge(i,j,k,c)/=coded_edge(gi,c,owner))bad=bad+1
              t=0
              do d=1,3
                if(d==c)cycle
                if(gi(d)>=0.and.gi(d)<nn(d)-1.and.mod(gi(d)+1,bn(d))==0)t=t+1
              end do
              if(t==2)corners=corners+1
            end do
          end do
        end do
      end do
      saved(g)%edge=fluxes(g)%edge
    end do
    call yy_sync_ct_edges(patch,fluxes)
    do n=1,igridstail
      g=igrids(n)
      if(any(fluxes(g)%edge/=saved(g)%edge))bad=bad+1
    end do
    call MPI_ALLREDUCE(corners,allcorners,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    if(allcorners==0)bad=bad+1

    ! Select only geometrical owners on angular interfaces, away from unresolved
    ! radial-wall edge stencils. Non-owner copies start deliberately corrupted.
    do pass=1,2
      t=0
      do n=1,igridstail
        g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
        do c=1,3
          fl=lo-1;fl(c)=lo(c)
          do k=fl(3),hi(3)
            do j=fl(2),hi(2)
              do i=fl(1),hi(1)
                idx=[i,j,k];gi=bc*bn+idx-lo
                if(.not.edge_selected(gi,c))cycle
                if(.not.yy_ct_edge_is_owner(g,idx,c))cycle
                t=t+1
                if(pass==2)receivers(:,t)=[g,idx,c]
              end do
            end do
          end do
        end do
      end do
      if(pass==1)allocate(receivers(5,t))
    end do
    call yy_prepare_ct_edges(cross,ps,receivers)
    do n=1,igridstail
      g=igrids(n)
      do c=1,3
        fl=lo-1;fl(c)=lo(c);offset=0.5d0;offset(c)=0.d0
        do k=fl(3),hi(3)
          do j=fl(2),hi(2)
            do i=fl(1),hi(1)
              idx=[i,j,k];pos=ps(g)%x(i,j,k,:)+offset*ps(g)%dx(i,j,k,:)
              call yy_basis(pos(2),pos(3),yy_patch,frame)
              value=0.37d0*ps(g)%dsC(i,j,k,c)*dot_product([0.2d0,-0.3d0,0.5d0],frame(:,c))
              if(.not.yy_ct_edge_is_owner(g,idx,c))value=value+7.d0
              fluxes(g)%edge(i,j,k,c)=value
            end do
          end do
        end do
      end do
    end do
    if(iprob==1)then
      ! Every second radial block has a non-owner edge on the physical theta
      ! boundary. Replacing the owner record must be rejected before exchange.
      do n=1,size(cross%receivers,2)
        g=cross%receivers(1,n);c=cross%receivers(5,n)
        if(node(pig1_,g)>1.and.c/=1)then
          cross%receivers(2,n)=lo(1)-1
          exit
        end if
      end do
    end if
    call yy_ordered_ct_edges(patch,cross,fluxes)
    if(iprob==1)call mpistop('test failed to select a non-owner receiver')
    do n=1,igridstail
      g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
      if(any(fluxes(g)%face/=123.d0))bad=bad+1
      do c=1,3
        fl=lo-1;fl(c)=lo(c);offset=0.5d0;offset(c)=0.d0
        do k=fl(3),hi(3)
          do j=fl(2),hi(2)
            do i=fl(1),hi(1)
              gi=bc*bn+[i,j,k]-lo
              if(.not.edge_selected(gi,c))cycle
              pos=ps(g)%x(i,j,k,:)+offset*ps(g)%dx(i,j,k,:)
              call yy_basis(pos(2),pos(3),yy_patch,frame)
              value=0.37d0*ps(g)%dsC(i,j,k,c)*dot_product([0.2d0,-0.3d0,0.5d0],frame(:,c))
              local(1)=max(local(1),abs(fluxes(g)%edge(i,j,k,c)-value))
            end do
          end do
        end do
      end do
    end do
    ! Directly compare every shared copy to its canonical source after overset.
    do n=1,size(patch%edge%send)
      g=patch%edge%request(1,n);idx=patch%edge%request(2:4,n)+lo;c=patch%edge%request(5,n)
      patch%edge%send(n)=fluxes(g)%edge(idx(1),idx(2),idx(3),c)
    end do
    call yy_patch_copy_exchange(patch%edge)
    do n=1,size(patch%edge%recv)
      g=patch%edge%receivers(1,n);idx=patch%edge%receivers(2:4,n)+lo;c=patch%edge%receivers(5,n)
      if(fluxes(g)%edge(idx(1),idx(2),idx(3),c)/=patch%edge%recv(n))bad=bad+1
    end do

    write(filename,'(a,i0,a)')'sync-rank-',yy_world_rank,'.txt';open(newunit=u,file=trim(filename),status='replace')
    do n=1,igridstail
      g=igrids(n)
      allocate(probe(g)%ws(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3))
      probe(g)%surfaceC=>ps(g)%surfaceC;probe(g)%dvolume=>ps(g)%dvolume
    end do
    mixed=0
    do field=0,1
      do n=1,igridstail
        g=igrids(n);probe(g)%ws=1.d20
        do c=1,3
          fl=lo;fl(c)=fl(c)-1
          do k=fl(3),hi(3)
            do j=fl(2),hi(2)
              do i=fl(1),hi(1)
                call face_value(g,[i,j,k],c,field,value)
                probe(g)%ws(i,j,k,c)=value
              end do
            end do
          end do
        end do
        base(g)%ws=probe(g)%ws
      end do
      call yy_sync_ct_faces(patch,probe,.false.)
      call yy_fill_ct_ghosts(ghost,probe,.true.)
      holes=0
      do n=1,igridstail
        g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
        do c=1,3
          fl=[ixGlo1,ixGlo2,ixGlo3];fl(c)=fl(c)-1
          do k=fl(3),ixGhi3
            do j=fl(2),ixGhi2
              do i=fl(1),ixGhi1
                gi=bc*bn+[i,j,k]-lo
                if(.not.radial_valid(gi,c))cycle
                if(abs(probe(g)%ws(i,j,k,c))>10.d0)holes=holes+1
              end do
            end do
          end do
        end do
      end do
      call MPI_ALLREDUCE(holes,allholes,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
      if(allholes==0)bad=bad+1 ! The mixed corners must actually exercise phase 3.
      call yy_sync_ct_faces(patch,probe,.true.)
      do n=1,igridstail
        g=igrids(n);manual(g)%ws=probe(g)%ws;probe(g)%ws=base(g)%ws
      end do
      do pass=1,3
        call yy_ordered_ct_ghosts(patch,ghost,probe)
        do n=1,igridstail
          g=igrids(n)
          if(any(probe(g)%ws/=manual(g)%ws))bad=bad+1
        end do
      end do
      do n=1,igridstail
        g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
        do c=1,3
          fl=[ixGlo1,ixGlo2,ixGlo3];fl(c)=fl(c)-1
          do k=fl(3),ixGhi3
            do j=fl(2),ixGhi2
              do i=fl(1),ixGhi1
                idx=[i,j,k];gi=bc*bn+idx-lo
                active=all(idx>=lo).and.all(idx<=hi)
                if(idx(c)==lo(c)-1)then
                  owner=idx;owner(c)=owner(c)+1
                  active=all(owner>=lo).and.all(owner<=hi)
                end if
                if(active.or..not.radial_valid(gi,c))then
                  if(probe(g)%ws(i,j,k,c)/=base(g)%ws(i,j,k,c))bad=bad+1
                end if
                if(.not.radial_valid(gi,c))cycle
                if(abs(probe(g)%ws(i,j,k,c))>10.d0)bad=bad+1
                write(u,'(6i8,es26.17)')field,yy_patch,gi,c,probe(g)%ws(i,j,k,c)
              end do
            end do
          end do
        end do
        block=>probe(g);old_stagger=stagger_grid;stagger_grid=.true.
        call get_divb(ps(g)%w,ixG^LL,ixG^LL,divb)
        stagger_grid=old_stagger;block=>ps(g)
        do k=ixGlo3,ixGhi3
          do j=ixGlo2,ixGhi2
            do i=ixGlo1,ixGhi1
              gi=bc*bn+[i,j,k]-lo
              if(gi(1)<0.or.gi(1)>=nn(1))cycle
              if(all(gi(2:3)>=0).and.all(gi(2:3)<nn(2:3)))cycle
              local(2)=max(local(2),abs(divb(i,j,k)))
              if(any([i,j,k]<lo).or.any([i,j,k]>hi))mixed=mixed+1
            end do
          end do
        end do
      end do
    end do
    close(u)
    call MPI_ALLREDUCE(local,maximum,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(bad,allbad,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(mixed,allmixed,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    if(yy_world_rank==0)write(*,'(a,4i10,2es26.17)')'CT_PATCH_SYNC ',allbad,allcorners,allholes,allmixed,maximum
    if(allbad/=0.or.maxval(maximum)>2.d-11)call mpistop('same-patch CT synchronization test failed')
    call sleep(1)
  contains
    double precision function coded_edge(gi,c,owner) result(v)
      integer,intent(in) :: gi(3),c,owner(3)
      v=yy_patch*1.d7+c*1.d6+(gi(1)+1)+101.d0*(gi(2)+1)+10007.d0*(gi(3)+1) &
        +1.d-4*(owner(1)+nb(1)*(owner(2)+nb(2)*owner(3)))
    end function
    logical function edge_selected(gi,c) result(ok)
      integer,intent(in) :: gi(3),c
      ok=gi(1)>=1.and.gi(1)<=nn(1)-2.and. &
        ((c/=2.and.(gi(2)==-1.or.gi(2)==nn(2)-1)).or. &
         (c/=3.and.(gi(3)==-1.or.gi(3)==nn(3)-1)))
    end function
    logical function radial_valid(gi,c) result(ok)
      integer,intent(in) :: gi(3),c
      integer :: minimum
      minimum=0;if(c==1)minimum=-1
      ok=gi(1)>=minimum.and.gi(1)<nn(1)
    end function
    subroutine face_value(g,idx,c,field,value)
      integer,intent(in) :: g,idx(3),c,field
      double precision,intent(out) :: value
      integer :: d,e,below(3)
      double precision :: upper,lower
      value=0.d0
      do d=1,3
        do e=1,3
          if(lvc(c,d,e)==0)cycle
          below=idx;below(d)=below(d)-1
          call potential_edge(g,idx,e,field,upper)
          call potential_edge(g,below,e,field,lower)
          value=value+lvc(c,d,e)*(upper-lower)
        end do
      end do
      value=value/ps(g)%surfaceC(idx(1),idx(2),idx(3),c)
    end subroutine
    subroutine potential_edge(g,idx,c,field,value)
      integer,intent(in) :: g,idx(3),c,field
      double precision,intent(out) :: value
      double precision,parameter :: nodes(4)=[-.861136311594052575d0,-.339981043584856265d0, &
        .339981043584856265d0,.861136311594052575d0]
      double precision,parameter :: weights(4)=[.347854845137453857d0,.652145154862546143d0, &
        .652145154862546143d0,.347854845137453857d0]
      integer :: above(3),q
      double precision :: off(3),pos(3),x(3),b(3,3),xyz(3),v(3),avec(3),jac,width
      above=max(idx,[ixGlo1,ixGlo2,ixGlo3]);off=.5d0;off(c)=0.d0
      off=off+dble(idx-above)
      pos=ps(g)%x(above(1),above(2),above(3),:)+off*ps(g)%dx(above(1),above(2),above(3),:)
      width=ps(g)%dx(above(1),above(2),above(3),c);value=0.d0
      do q=1,4
        x=pos;x(c)=x(c)+.5d0*width*nodes(q)
        call yy_basis(x(2),x(3),yy_patch,b);xyz=x(1)*b(:,1)
        v=[.2d0,-.3d0,.5d0];if(field==1)v=[xyz(2),xyz(3),xyz(1)]
        avec=[v(2)*xyz(3)-v(3)*xyz(2),v(3)*xyz(1)-v(1)*xyz(3),v(1)*xyz(2)-v(2)*xyz(1)]
        avec=avec/2.d0;if(field==1)avec=avec*(2.d0/3.d0)
        jac=1.d0;if(c==2)jac=x(1)
        if(c==3)jac=x(1)*sin(x(2))
        value=value+.5d0*width*weights(q)*jac*dot_product(avec,b(:,c))
      end do
    end subroutine
  end subroutine
end module
