program check_ct_coupling
  use mpi
  use mod_yinyang_staggered, only: yy_staggered_basis,yy_staggered_offset
  use mod_yinyang_ct_coupling
  implicit none
  type(yy_ct_edge_plan) :: plan
  integer :: ierr,rank,np,patch,scale,kind,mode,nn(3),bn(3),nb(3),n,p,c,i,j,k,h,seq,nt,pass,q,g,caseid,bad,badall,total
  integer,allocatable :: directory(:,:),directions(:),ids(:)
  double precision,allocatable :: radial(:),points(:,:),widths(:),donor(:),answer(:),saved(:)
  double precision :: pi,angular(2,2),step(2),off(3),coord(3),frame(3,3),xyz(3),vec(3),width,length,dt,ref,err,sq,emax,rms
  character(len=128) :: arg,filename
  integer :: u
  call MPI_INIT(ierr);call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr);call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  call get_command_argument(1,arg);read(arg,*)scale
  call get_command_argument(2,arg);read(arg,*)kind
  call get_command_argument(3,arg);read(arg,*)mode
  patch=rank/(np/2);pi=acos(-1.d0);nn=[8,16,48]*scale;bn=4;nb=nn/bn
  allocate(radial(0:nn(1)),directory(6,2*product(nb)))
  do i=0,nn(1)
    radial(i)=1.d0+dble(i)/nn(1)
    if(mode==1)radial(i)=2.d0**(dble(i)/nn(1))
  end do
  angular(:,1)=[0.21875d0*pi,-0.78125d0*pi];angular(:,2)=[0.78125d0*pi,0.78125d0*pi]
  step=(angular(:,2)-angular(:,1))/nn(2:3)
  n=0
  do p=0,1
    do k=0,nb(3)-1
      do j=0,nb(2)-1
        do i=0,nb(1)-1
          n=n+1;directory(:,n)=[i,j,k,n,p*(np/2)+mod(n-1,np/2),p]
        end do
      end do
    end do
  end do
  do pass=1,2
    n=0;seq=0
    do c=1,3
      call yy_staggered_offset(2,c,off)
      do h=1,3
        i=h*nn(1)/4
        do k=-1,nn(3)
          do j=-1,nn(2)
            if(j>=0.and.j<nn(2).and.k>=0.and.k<nn(3)) cycle
            if(j>=0.and.j<nn(2).and.mod(j,8)/=0) cycle
            if(k>=0.and.k<nn(3).and.mod(k,8)/=0) cycle
            seq=seq+1
            if(mod(seq-1,np/2)/=mod(rank,np/2)) cycle
            n=n+1
            if(pass==2) then
              coord=[radial(i+1),angular(:,1)+([j,k]+off(2:3))*step]
              if(c==1)coord(1)=0.5d0*(radial(i)+radial(i+1))
              points(:,n)=coord;directions(n)=c;ids(n)=seq
              widths(n)=radial(i+1)-radial(i)
              if(c>1)widths(n)=step(c-1)
            end if
          end do
        end do
      end do
    end do
    if(pass==1)then
      nt=n+1;allocate(points(3,nt),widths(nt),directions(nt),ids(nt),answer(nt),saved(nt))
    end if
  end do
  points(:,nt)=[0.5d0,pi/2,0.d0];widths(nt)=0.01d0;directions(nt)=1;ids(nt)=-1
  call yy_ct_edge_build(plan,MPI_COMM_WORLD,patch,bn,nn,radial,angular,directory,points,widths,directions)
  bad=0
  if(plan%covered(nt).or.any(.not.plan%covered(:nt-1)))bad=bad+1
  allocate(donor(size(plan%samples%request,2)));dt=0.37d0
  write(filename,'(a,i0,a)')'rank-',rank,'.txt';open(newunit=u,file=trim(filename),status='replace')
  do caseid=0,1
    do n=1,size(donor)
      g=plan%samples%request(1,n);c=plan%samples%request(5,n)
      i=directory(1,g)*bn(1)+plan%samples%request(2,n)
      j=directory(2,g)*bn(2)+plan%samples%request(3,n)
      k=directory(3,g)*bn(3)+plan%samples%request(4,n)
      call yy_staggered_offset(2,c,off)
      coord=[radial(i+1),angular(:,1)+([j,k]+off(2:3))*step]
      width=radial(i+1)-radial(i)
      if(c==1)coord(1)=0.5d0*(radial(i)+radial(i+1))
      if(c>1)width=step(c-1)
      length=edge_length(coord,width,c)
      if(abs(length-plan%donor_length(n))>1.d-13)bad=bad+1
      call yy_staggered_basis(coord(2),coord(3),patch,frame)
      xyz=coord(1)*frame(:,1);call analytic(caseid,xyz,vec)
      donor(n)=dt*length*dot_product(vec,frame(:,c))
      if(kind==1)donor(n)=dt*integral(caseid,patch,c,coord,width)
    end do
    do pass=1,3
      call yy_ct_edge_exchange(plan,donor,answer)
      if(pass==1)saved=answer
      if(any(answer/=saved))bad=bad+1
    end do
    err=0.d0;sq=0.d0
    do n=1,nt-1
      ref=dt*integral(caseid,patch,directions(n),points(:,n),widths(n))
      length=edge_length(points(:,n),widths(n),directions(n))
      err=max(err,abs(answer(n)-ref)/(dt*length));sq=sq+((answer(n)-ref)/(dt*length))**2
      write(u,'(3i9,es26.17)')caseid,patch,ids(n),answer(n)
    end do
    if(answer(nt)/=0.d0)bad=bad+1
    call MPI_ALLREDUCE(err,emax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(sq,rms,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    n=nt-1;call MPI_ALLREDUCE(n,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(rank==0)write(*,'(a,5i8,2es26.17)')'RESULT ',scale,kind,mode,caseid,total,emax,sqrt(rms/total)
    if(kind==0.and.caseid==0.and.emax>1.d-7)bad=bad+1
    donor=1.7d0*donor;call yy_ct_edge_exchange(plan,donor,answer)
    if(maxval(abs(answer-1.7d0*saved))>1.d-13)bad=bad+1
    donor=0.d0;call yy_ct_edge_exchange(plan,donor,answer)
    if(any(answer/=0.d0))bad=bad+1
  end do
  close(u)
  ! Empty receiver ranks still supply donor data collectively.
  call yy_ct_edge_clear(plan)
  n=0;if(rank==0)n=1
  call yy_ct_edge_build(plan,MPI_COMM_WORLD,patch,bn,nn,radial,angular,directory, &
    points(:,:n),widths(:n),directions(:n))
  deallocate(donor,answer);allocate(donor(size(plan%donor_length)),answer(n))
  donor=0.d0;call yy_ct_edge_exchange(plan,donor,answer)
  if(any(answer/=0.d0))bad=bad+1
  call yy_ct_edge_clear(plan)
  call ghost_test(bad)
  call MPI_ALLREDUCE(bad,badall,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(rank==0)write(*,'(a,i8)')'CHECK ',badall
  if(badall/=0)call MPI_ABORT(MPI_COMM_WORLD,41,ierr)
  call MPI_FINALIZE(ierr)
contains
  pure subroutine analytic(which,x,v)
    integer,intent(in)::which
    double precision,intent(in)::x(3)
    double precision,intent(out)::v(3)
    v=[0.2d0,-0.3d0,0.5d0]
    if(which==1)v=[sin(x(2)),sin(x(3)),sin(x(1))]
  end subroutine
  pure double precision function edge_length(pos,width,c)result(v)
    double precision,intent(in)::pos(3),width
    integer,intent(in)::c
    v=width
    if(c==2)v=width*pos(1)
    if(c==3)v=width*pos(1)*sin(pos(2))
  end function
  ! Independent composite Simpson reference (128 panels), including rotation.
  double precision function integral(which,p,c,pos,width)result(v)
    integer,intent(in)::which,p,c
    double precision,intent(in)::pos(3),width
    double precision::x(3),e(3),r,t,f,b(3),jac,weight
    integer::a
    v=0.d0
    do a=0,128
      r=pos(1);t=pos(2);f=pos(3)
      if(c==1)r=r+width*(dble(a)/128-0.5d0)
      if(c==2)t=t+width*(dble(a)/128-0.5d0)
      if(c==3)f=f+width*(dble(a)/128-0.5d0)
      x=[r*sin(t)*cos(f),r*sin(t)*sin(f),r*cos(t)]
      select case(c)
      case(1);e=x/r;jac=1.d0
      case(2);e=[cos(t)*cos(f),cos(t)*sin(f),-sin(t)];jac=r
      case(3);e=[-sin(f),cos(f),0.d0];jac=r*sin(t)
      end select
      if(p==1)then
        x=[-x(1),x(3),x(2)];e=[-e(1),e(3),e(2)]
      end if
      call analytic(which,x,b)
      weight=2.d0
      if(mod(a,2)==1)weight=4.d0
      if(a==0.or.a==128)weight=1.d0
      v=v+weight*dot_product(e,b)*jac
    end do
    v=v*width/(3*128)
  end function
  subroutine ghost_test(bad)
    integer,intent(inout)::bad
    double precision::area(0:8,0:8,0:8,3),bf(0:8,0:8,0:8,3),old(0:8,0:8,0:8,3)
    double precision::after(0:8,0:8,0:8,3),rr(0:8),tt(0:8),pp(0:8),balance,peak,vol
    logical::sides(2,2),selected,allowed(0:8,0:8,0:8,3)
    integer::i,j,k,c,mask,status
    do i=0,8
      rr(i)=1.d0+0.1d0*i
      if(mode==1)rr(i)=1.08d0**i
      tt(i)=0.7d0+0.06d0*i;pp(i)=-0.3d0+0.05d0*i
    end do
    area=1.d0
    do k=1,8
      do j=1,8
        do i=1,8
          area(i,j,k,1)=rr(i)**2*(cos(tt(j-1))-cos(tt(j)))*(pp(k)-pp(k-1))
          area(i,j,k,2)=0.5d0*(rr(i)**2-rr(i-1)**2)*sin(tt(j))*(pp(k)-pp(k-1))
          area(i,j,k,3)=0.5d0*(rr(i)**2-rr(i-1)**2)*(tt(j)-tt(j-1))
        end do
      end do
    end do
    area(:,0,:,2)=area(:,1,:,2)*sin(tt(0))/sin(tt(1))
    area(:,:,0,3)=area(:,:,1,3)
    do c=1,3
      do k=0,8
        do j=0,8
          do i=0,8
            old(i,j,k,c)=sin(0.2d0*i+0.3d0*j-0.4d0*k+c)
          end do
        end do
      end do
    end do
    peak=0.d0
    ! All 16 combinations, including corners and interior blocks (no sides).
    do mask=0,15
      sides=reshape([(btest(mask,c),c=0,3)],[2,2]);bf=old;allowed=.false.
      if(sides(1,1))allowed(3:6,0:1,3:6,2)=.true.
      if(sides(2,1))allowed(3:6,7:8,3:6,2)=.true.
      j=3;k=6
      if(sides(1,1))j=1
      if(sides(2,1))k=8
      if(sides(1,2))allowed(3:6,j:k,0:1,3)=.true.
      if(sides(2,2))allowed(3:6,j:k,7:8,3)=.true.
      call yy_ct_close_ghosts([0,0,0],[3,3,3],[6,6,6],[1,1,1],[8,8,8],sides,bf,area,status)
      if(status/=0)bad=bad+1
      if(any(bf/=old.and..not.allowed))bad=bad+1
      do k=1,8
        do j=1,8
          selected=((j<3.and.sides(1,1)).or.(j>6.and.sides(2,1))).and.k>=3.and.k<=6
          if((k<3.and.sides(1,2)).or.(k>6.and.sides(2,2))) &
            selected=(j>=3.and.j<=6).or.(j<3.and.sides(1,1)).or.(j>6.and.sides(2,1))
          if(.not.selected)cycle
          do i=3,6
            balance=area(i,j,k,1)*bf(i,j,k,1)-area(i-1,j,k,1)*bf(i-1,j,k,1) &
              +area(i,j,k,2)*bf(i,j,k,2)-area(i,j-1,k,2)*bf(i,j-1,k,2) &
              +area(i,j,k,3)*bf(i,j,k,3)-area(i,j,k-1,3)*bf(i,j,k-1,3)
            vol=(rr(i)**3-rr(i-1)**3)/3*(cos(tt(j-1))-cos(tt(j)))*(pp(k)-pp(k-1))
            peak=max(peak,abs(balance)/vol)
          end do
        end do
      end do
      after=bf
      call yy_ct_close_ghosts([0,0,0],[3,3,3],[6,6,6],[1,1,1],[8,8,8],sides,bf,area,status)
      if(any(bf/=after))bad=bad+1
    end do
    if(peak>1.d-11)bad=bad+1
    ! Invalid geometry must fail before modifying a single face.
    bf=old;area(3,1,3,2)=0.d0;sides=.true.
    call yy_ct_close_ghosts([0,0,0],[3,3,3],[6,6,6],[1,1,1],[8,8,8],sides,bf,area,status)
    if(status==0.or.any(bf/=old))bad=bad+1
    if(rank==0)write(*,'(a,es26.17)')'GHOST_DIVB ',peak
  end subroutine
end program
