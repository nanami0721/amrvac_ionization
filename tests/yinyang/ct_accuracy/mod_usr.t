! Isolate CT transfer errors on the actual mesh without modifying evolved states.
module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved
  use mod_yinyang, only: yy_basis,yy_ct_patch_plan,yy_prepare_ct_patch,yy_ct_stage_plan, &
    yy_prepare_ct_edges,yy_ordered_ct_edges,yy_ct_edge_is_owner,yy_ct_ghost_plan, &
    yy_prepare_ct_ghosts,yy_fill_ct_ghosts
  use mod_yinyang_ct_coupling, only: yy_ct_close_ghosts
  use mod_functions_bfield, only: get_divb
  implicit none
  logical  ::  checked=.false.
  integer  ::  field_id=0,quadrature_order=8
  double precision  ::  gx(16),gw(16),velocity(3),normal(3),origin(3),radius=.65d0,amplitude
contains
  subroutine usr_init()
    integer  ::  n
    namelist /accuracy_test/ quadrature_order
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,accuracy_test,end=10)
10    close(unitpar)
    end do
    if(quadrature_order/=8.and.quadrature_order/=16)call mpistop('quadrature order must be 8 or 16')
    call gauss_rule(quadrature_order)
    call set_coordinate_system('yinyang');call mhd_activate()
    usr_init_one_grid=>initial;usr_process_global=>diagnostic
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    if(allocated(block%ws))block%ws=0.d0
    w(ixO^S,:)=0.d0;w(ixO^S,rho_)=1.d0;w(ixO^S,p_)=1.d0
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine
  pure function cross(a,b) result(c)
    double precision,intent(in) :: a(3),b(3)
    double precision :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function
  pure subroutine vectors(pos,avec,bvec)
    double precision,intent(in) :: pos(3)
    double precision,intent(out) :: avec(3),bvec(3)
    double precision :: q(3),s,a0
    select case(field_id)
    case(0)
      bvec=[.2d0,-.1d0,.3d0];avec=cross(bvec,pos)/2.d0
    case(1)
      bvec=.2d0*[pos(2),pos(3),pos(1)];avec=cross(bvec,pos)/3.d0
    case default
      q=pos-origin;s=sum(q*q)/radius**2;avec=0.d0;bvec=0.d0
      if(s<1.d0)then
        a0=amplitude*radius/((8.d0/sqrt(7.d0))*(6.d0/7.d0)**3)
        avec=a0*(1.d0-s)**4*normal
        bvec=-8.d0*a0/radius**2*(1.d0-s)**3*cross(q,normal)
      end if
    end select
  end subroutine
  pure subroutine location(g,idx,c,coord,width)
    integer,intent(in) :: g,idx(3),c
    double precision,intent(out) :: coord(3),width(3)
    integer :: gi(3)
    double precision :: fraction(3)
    gi=([node(pig1_,g),node(pig2_,g),node(pig3_,g)]-1)*[block_nx1,block_nx2,block_nx3]+ &
      idx-[ixMlo1,ixMlo2,ixMlo3]
    width=([xprobmax1,xprobmax2,xprobmax3]-[xprobmin1,xprobmin2,xprobmin3])/ &
      [domain_nx1,domain_nx2,domain_nx3]
    fraction=1.d0;fraction(c)=.5d0
    coord=[xprobmin1,xprobmin2,xprobmin3]+(gi+fraction)*width
  end subroutine
  ! mode 0 is midpoint E; 1 is integrated E; 2 is integrated vector potential.
  pure subroutine edge_integral(g,idx,c,mode,value)
    integer,intent(in) :: g,idx(3),c,mode
    double precision,intent(out) :: value
    double precision :: coord(3),x(3),width(3),frame(3,3),a(3),b(3),vec(3),jac,w,point
    integer :: q,nq
    call location(g,idx,c,coord,width)
    nq=quadrature_order;if(mode==0)nq=1
    value=0.d0
    do q=1,nq
      point=0.d0;w=2.d0
      if(mode/=0)then
        point=gx(q);w=gw(q)
      end if
      x=coord;x(c)=x(c)+.5d0*width(c)*point
      call yy_basis(x(2),x(3),yy_patch,frame)
      call vectors(x(1)*frame(:,1),a,b)
      vec=-cross(velocity,b);if(mode==2)vec=a
      jac=1.d0
      if(c==2)jac=x(1)
      if(c==3)jac=x(1)*sin(x(2))
      value=value+.5d0*width(c)*w*jac*dot_product(vec,frame(:,c))
    end do
  end subroutine
  pure function face_reference(g,idx,c) result(value)
    integer,intent(in) :: g,idx(3),c
    integer :: d,e,below(3)
    double precision :: value,upper_edge,lower_edge
    value=0.d0
    do d=1,3
      do e=1,3
        if(lvc(c,d,e)==0)cycle
        below=idx;below(d)=below(d)-1
        call edge_integral(g,idx,e,2,upper_edge)
        call edge_integral(g,below,e,2,lower_edge)
        value=value+lvc(c,d,e)*(upper_edge-lower_edge)
      end do
    end do
    value=value/ps(g)%surfaceC(idx(1),idx(2),idx(3),c)
  end function
  subroutine diagnostic(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    type(yy_ct_patch_plan) :: patch
    type(yy_ct_stage_plan) :: edges,legacy
    type(yy_ct_ghost_plan) :: ghost
    type(ct_stage_flux),allocatable :: mapped(:),midpoint(:),exactedge(:)
    type(state),allocatable,target,save :: reference(:),probe(:)
    integer,allocatable :: receivers(:,:)
    integer :: n,g,c,i,j,k,t,pass,mode,stage,d,e,status,u,lo(3),hi(3),fl(3),fh(3),idx(3),gi(3),bc(3),bn(3),nn(3),nb(3),below(3)
    double precision :: v(14),total(14),mx(6),maximum(6),ref,raw,err,midcurl,newcurl,exactcurl,area,ds,divb(ixG^T),phi,tm
    logical :: sides(2,2),seam
    character(len=80) :: filename
    if(checked)return
    checked=.true.;lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3]
    bn=[block_nx1,block_nx2,block_nx3];nn=[domain_nx1,domain_nx2,domain_nx3];nb=nn/bn
    allocate(mapped(max_blocks),midpoint(max_blocks),exactedge(max_blocks),reference(max_blocks),probe(max_blocks))
    do n=1,igridstail
      g=igrids(n)
      allocate(mapped(g)%edge(ixG^T,3),midpoint(g)%edge(ixG^T,3),exactedge(g)%edge(ixG^T,3))
      allocate(reference(g)%ws(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3))
      allocate(probe(g)%ws(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3))
      reference(g)%surfaceC=>ps(g)%surfaceC;reference(g)%dvolume=>ps(g)%dvolume
      probe(g)%surfaceC=>ps(g)%surfaceC;probe(g)%dvolume=>ps(g)%dvolume
    end do
    call yy_prepare_ct_patch(patch,ps);call yy_prepare_ct_ghosts(ghost,ps,face_moments=.true.)
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
                if(.not.((c/=2.and.(gi(2)==-1.or.gi(2)==nn(2)-1)).or. &
                  (c/=3.and.(gi(3)==-1.or.gi(3)==nn(3)-1))))cycle
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
    call yy_prepare_ct_edges(edges,ps,receivers,radial_extrapolate=.true.)
    call yy_prepare_ct_edges(legacy,ps,receivers,radial_extrapolate=.true.,native_midpoint=.false.)
    call sleep(1)
    write(filename,'(a,i0,a)')'accuracy-rank-',yy_world_rank,'.txt'
    open(newunit=u,file=trim(filename),status='replace')
    do field_id=0,3
      velocity=[.4d0,-.2d0,.3d0];amplitude=1.d0
      if(field_id>=2)then
        phi=146.25d0*acos(-1.d0)/180.d0;normal=[cos(phi),sin(phi),0.d0]
        velocity=[-sin(phi),cos(phi),0.d0];tm=.404d0
        if(field_id==3)then
          velocity=-velocity;tm=.522d0
        end if
        origin=1.5d0*normal+(tm-.3d0)*velocity;amplitude=1.d-3
      end if
      do n=1,igridstail
        g=igrids(n);midpoint(g)%edge=0.d0;exactedge(g)%edge=0.d0;reference(g)%ws=0.d0
        do c=1,3
          fl=lo-1;fl(c)=lo(c)
          do k=fl(3),hi(3)
            do j=fl(2),hi(2)
              do i=fl(1),hi(1)
                call edge_integral(g,[i,j,k],c,0,midpoint(g)%edge(i,j,k,c))
                call edge_integral(g,[i,j,k],c,1,exactedge(g)%edge(i,j,k,c))
              end do
            end do
          end do
          fl=[ixMlo1,ixGlo2,ixGlo3];fh=[ixMhi1,ixGhi2,ixGhi3];fl(c)=fl(c)-1
          do k=fl(3),fh(3)
            do j=fl(2),fh(2)
              do i=fl(1),fh(1)
                reference(g)%ws(i,j,k,c)=face_reference(g,[i,j,k],c)
              end do
            end do
          end do
        end do
      end do
      do mode=0,3
        do n=1,igridstail
          g=igrids(n);mapped(g)%edge=midpoint(g)%edge
          if(mode==1)mapped(g)%edge=exactedge(g)%edge
        end do
        if(mode==3)then
          call yy_ordered_ct_edges(patch,edges,mapped)
        else
          call yy_ordered_ct_edges(patch,legacy,mapped)
        end if
        ! Oracle control only: express receiver integrals in the native midpoint
        ! convention. This uses analytic data and is NOT a production algorithm.
        if(mode==2)then
          do n=1,igridstail
            g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
            do c=1,3
              fl=lo-1;fl(c)=lo(c)
              do k=fl(3),hi(3)
                do j=fl(2),hi(2)
                  do i=fl(1),hi(1)
                    gi=bc*bn+[i,j,k]-lo
                    if(.not.((c/=2.and.(gi(2)==-1.or.gi(2)==nn(2)-1)).or. &
                      (c/=3.and.(gi(3)==-1.or.gi(3)==nn(3)-1))))cycle
                    mapped(g)%edge(i,j,k,c)=mapped(g)%edge(i,j,k,c) &
                      +midpoint(g)%edge(i,j,k,c)-exactedge(g)%edge(i,j,k,c)
                  end do
                end do
              end do
            end do
          end do
        end if
        v=0.d0;mx=0.d0
        do n=1,size(receivers,2)
          g=receivers(1,n);idx=receivers(2:4,n);c=receivers(5,n)
          ds=ps(g)%dsC(idx(1),idx(2),idx(3),c)
          ref=exactedge(g)%edge(idx(1),idx(2),idx(3),c)/ds
          raw=midpoint(g)%edge(idx(1),idx(2),idx(3),c)/ds
          err=mapped(g)%edge(idx(1),idx(2),idx(3),c)/ds-ref
          v(1:4)=v(1:4)+[err**2,ref**2,(raw-ref)**2,1.d0]
          mx(1)=max(mx(1),abs(err)/amplitude)
          gi=([node(pig1_,g),node(pig2_,g),node(pig3_,g)]-1)*bn+idx-lo
          write(u,'(a,7i7,3es27.17e3)')'E',field_id,mode,yy_patch,gi,c,ref,raw,mapped(g)%edge(idx(1),idx(2),idx(3),c)/ds
        end do
        do n=1,igridstail
          g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
          do c=1,3
            fl=lo;fl(c)=lo(c)-1
            do k=fl(3),hi(3)
              do j=fl(2),hi(2)
                do i=fl(1),hi(1)
                  idx=[i,j,k];gi=bc*bn+idx-lo;midcurl=0.d0;newcurl=0.d0;exactcurl=0.d0
                  do d=1,3
                    do e=1,3
                      if(lvc(c,d,e)==0)cycle
                      below=idx;below(d)=below(d)-1
                      midcurl=midcurl+lvc(c,d,e)*(midpoint(g)%edge(i,j,k,e)-midpoint(g)%edge(below(1),below(2),below(3),e))
                      newcurl=newcurl+lvc(c,d,e)*(mapped(g)%edge(i,j,k,e)-mapped(g)%edge(below(1),below(2),below(3),e))
                      exactcurl=exactcurl+lvc(c,d,e)*(exactedge(g)%edge(i,j,k,e)-exactedge(g)%edge(below(1),below(2),below(3),e))
                    end do
                  end do
                  area=ps(g)%surfaceC(i,j,k,c)
                  midcurl=midcurl/area;newcurl=newcurl/area;exactcurl=exactcurl/area
                  v(5:8)=v(5:8)+[(newcurl-exactcurl)**2,exactcurl**2,(midcurl-exactcurl)**2,(newcurl-midcurl)**2]
                  mx(2)=max(mx(2),abs(newcurl-midcurl)/amplitude)
                  seam=any(gi(2:3)<=0).or.any(gi(2:3)>=nn(2:3)-1)
                  if(.not.seam)mx(3)=max(mx(3),abs(newcurl-midcurl)/amplitude)
                  if(seam)v(9:12)=v(9:12)+[(newcurl-exactcurl)**2,exactcurl**2,(midcurl-exactcurl)**2,1.d0]
                end do
              end do
            end do
          end do
        end do
        call MPI_ALLREDUCE(v,total,14,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
        call MPI_ALLREDUCE(mx,maximum,6,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
        if(yy_world_rank==0)write(*,'(a,2i4,20es27.17e3)')'EDGE_ACCURACY ',field_id,mode,total,maximum
      end do
      ! Stage 0 is compatible reference+closure, 1 interpolation, 2 interpolation+closure.
      do stage=0,2
        do n=1,igridstail
          g=igrids(n);probe(g)%ws=reference(g)%ws
        end do
        if(stage>0)call yy_fill_ct_ghosts(ghost,probe,.false.)
        if(stage/=1)then
          do n=1,igridstail
            g=igrids(n);bc=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
            sides(1,:)=(bc(2:3)==1);sides(2,:)=(bc(2:3)==nb(2:3))
            call yy_ct_close_ghosts([ixGlo1-1,ixGlo2-1,ixGlo3-1],lo,hi,[ixGlo1,ixGlo2,ixGlo3], &
              [ixGhi1,ixGhi2,ixGhi3],sides,probe(g)%ws,ps(g)%surfaceC,status)
            if(status/=0)call mpistop('analytic ghost closure rejected state')
          end do
        end if
        v=0.d0;mx=0.d0
        do n=1,size(ghost%values)
          g=ghost%receivers(1,n);idx=ghost%receivers(2:4,n);c=ghost%receivers(5,n)
          ref=reference(g)%ws(idx(1),idx(2),idx(3),c);raw=probe(g)%ws(idx(1),idx(2),idx(3),c);err=raw-ref
          v(1:3)=v(1:3)+[err**2,ref**2,1.d0];mx(1)=max(mx(1),abs(err)/amplitude)
          v(3*c+1:3*c+3)=v(3*c+1:3*c+3)+[err**2,ref**2,1.d0]
          gi=([node(pig1_,g),node(pig2_,g),node(pig3_,g)]-1)*bn+idx-lo
          write(u,'(a,7i7,2es27.17e3)')'B',field_id,stage,yy_patch,gi,c,ref,raw
        end do
        do n=1,igridstail
          g=igrids(n)
          do c=1,3
            fl=lo;fl(c)=lo(c)-1
            if(any(probe(g)%ws(fl(1):hi(1),fl(2):hi(2),fl(3):hi(3),c)/= &
              reference(g)%ws(fl(1):hi(1),fl(2):hi(2),fl(3):hi(3),c)))call mpistop('changed active face')
          end do
          block=>probe(g)
          call get_divb(ps(g)%w,ixG^LL,ixG^LL,divb)
          bc=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
          do k=ixGlo3,ixGhi3
            do j=ixGlo2,ixGhi2
              seam=((j<lo(2).and.bc(2)==1).or.(j>hi(2).and.bc(2)==nb(2))).and.k>=lo(3).and.k<=hi(3)
              if((k<lo(3).and.bc(3)==1).or.(k>hi(3).and.bc(3)==nb(3))) &
                seam=(j>=lo(2).and.j<=hi(2)).or.(j<lo(2).and.bc(2)==1).or.(j>hi(2).and.bc(2)==nb(2))
              if(seam)mx(2)=max(mx(2),maxval(abs(divb(lo(1):hi(1),j,k)))/amplitude)
            end do
          end do
          block=>ps(g)
        end do
        call MPI_ALLREDUCE(v,total,14,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
        call MPI_ALLREDUCE(mx,maximum,6,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
        if(yy_world_rank==0)write(*,'(a,2i4,20es27.17e3)')'GHOST_ACCURACY ',field_id,stage,total,maximum
        if(stage/=1.and.maximum(2)>1.d-10)call mpistop('constrained analytic ghost divergence')
        if(stage==0.and.maximum(1)>1.d-11)call mpistop('closure altered compatible reference')
      end do
    end do
    close(u)
  end subroutine
  subroutine gauss_rule(n)
    integer,intent(in) :: n
    integer :: i,j,iteration
    double precision :: z,old,p1,p2,p3,deriv
    do i=1,(n+1)/2
      z=cos(acos(-1.d0)*(i-.25d0)/(n+.5d0))
      do iteration=1,100
        p1=1.d0;p2=0.d0
        do j=1,n
          p3=p2;p2=p1;p1=((2*j-1)*z*p2-(j-1)*p3)/j
        end do
        deriv=n*(z*p1-p2)/(z*z-1.d0);old=z;z=old-p1/deriv
        if(abs(z-old)<2.d-15)exit
      end do
      gx(i)=-z;gx(n+1-i)=z;gw(i)=2.d0/((1.d0-z*z)*deriv**2);gw(n+1-i)=gw(i)
    end do
  end subroutine
end module
