program check_schedule
  use mpi
  use mod_yinyang_staggered
  implicit none
  type(yy_staggered_table) :: table
  integer :: ierr,rank,np,patch,scale,family,mode,nn(3),bn(3),nb(3),nblocks,n,p,a,b,c,i,j,k,h,side
  integer :: g,owner,nt,seq,target,u,caseid,repeat,covered,total,bad,badall
  integer,allocatable :: directory(:,:),indices(:)
  double precision,allocatable :: radial(:),points(:,:),samples(:),value(:,:),saved(:,:)
  double precision :: angular(2,2),step(2),off(3),coord(3),frame(3,3),xyz(3),vector(3),err,peak,geom,geomall
  double precision :: sumerr,globalerr,exact(3),r,pi
  character(len=128) :: arg,filename
  call MPI_INIT(ierr);call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr);call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
  call get_command_argument(1,arg);read(arg,*) scale
  call get_command_argument(2,arg);read(arg,*) family
  call get_command_argument(3,arg);read(arg,*) mode
  if(mod(np,2)/=0) stop 1
  patch=rank/(np/2);nn=[8*scale,16*scale,48*scale];bn=4;nb=nn/bn;pi=acos(-1.d0)
  angular(:,1)=[0.21875d0*pi,-0.78125d0*pi];angular(:,2)=[0.78125d0*pi,0.78125d0*pi]
  step=(angular(:,2)-angular(:,1))/nn(2:3)
  allocate(radial(0:nn(1)))
  do i=0,nn(1)
    radial(i)=1.d0+dble(i)/nn(1)
    if(mode==1) radial(i)=2.d0**(dble(i)/nn(1))
  end do
  nblocks=2*product(nb);allocate(directory(6,nblocks));n=0
  do p=0,1
    do k=0,nb(3)-1
      do j=0,nb(2)-1
        do i=0,nb(1)-1
          n=n+1;owner=p*(np/2)+mod(n-1,np/2)
          directory(:,n)=[i,j,k,n,owner,p]
        end do
      end do
    end do
  end do
  ! Same physical receivers for any MPI size; only ownership changes. Sample
  ! all three face/edge orientations, both seams and corners, three shell radii.
  nt=0;seq=0
  do repeat=1,2
    target=0;seq=0
    do c=1,3
      call yy_staggered_offset(family,c,off)
      do h=1,3
        i=h*nn(1)/4
        do k=-1,nn(3)
          do j=-1,nn(2)
            if(j>=0.and.j<nn(2).and.k>=0.and.k<nn(3)) cycle
            seq=seq+1
            if(mod(seq-1,np/2)/=mod(rank,np/2)) cycle
            target=target+1
            if(repeat==2) then
              r=radial(i+1)
              if(off(1)==0.5d0) r=0.5d0*(radial(i)+radial(i+1))
              points(:,target)=[r,angular(:,1)+([j,k]+off(2:3))*step]
              indices(target)=seq
            end if
          end do
        end do
      end do
    end do
    if(repeat==1) then
      nt=target+2;allocate(points(3,nt),indices(nt))
    end if
  end do
  ! Deliberately uncovered radial point, must not trigger invalid donors.
  points(:,nt-1)=[0.5d0,pi/2,0.d0];indices(nt-1)=-1
  points(:,nt)=[1.5d0,pi/2,0.d0];indices(nt)=-2
  call yy_staggered_build(table,MPI_COMM_WORLD,patch,family,bn,nn,radial,angular,directory,points)
  allocate(samples(size(table%request,2)),value(3,nt),saved(3,nt))
  geom=0.d0;bad=0
  if(any(table%covered(nt-1:nt))) bad=bad+1
  if(any(.not.table%covered(:nt-2))) bad=bad+1
  do n=1,size(samples)
    g=table%request(1,n);c=table%request(5,n)
    if(directory(5,g)/=rank.or.directory(6,g)/=patch) bad=bad+1
    call yy_staggered_offset(family,c,off)
    i=directory(1,g)*bn(1)+table%request(2,n)
    j=directory(2,g)*bn(2)+table%request(3,n)
    k=directory(3,g)*bn(3)+table%request(4,n)
    r=radial(i+1)
    if(off(1)==0.5d0) r=0.5d0*(radial(i)+radial(i+1))
    coord=[r,angular(:,1)+([j,k]+off(2:3))*step]
    call reference_frame(coord(2),coord(3),patch,frame)
    xyz=r*frame(:,1)
    geom=max(geom,maxval(abs(xyz-table%donor_position(:,n))),maxval(abs(frame(:,c)-table%donor_basis(:,n))))
  end do
  write(filename,'(a,i0,a)') 'rank-',rank,'.txt'
  open(newunit=u,file=trim(filename),status='replace')
  do caseid=0,1
    do n=1,size(samples)
      ! Independently decode the returned block/index requests, not cached xyz.
      g=table%request(1,n);c=table%request(5,n)
      call yy_staggered_offset(family,c,off)
      i=directory(1,g)*bn(1)+table%request(2,n)
      j=directory(2,g)*bn(2)+table%request(3,n)
      k=directory(3,g)*bn(3)+table%request(4,n)
      r=radial(i+1)
      if(off(1)==0.5d0) r=0.5d0*(radial(i)+radial(i+1))
      coord=[r,angular(:,1)+([j,k]+off(2:3))*step]
      call reference_frame(coord(2),coord(3),patch,frame)
      xyz=r*frame(:,1);call analytic(caseid,xyz,vector)
      samples(n)=dot_product(vector,frame(:,c))
    end do
    do repeat=1,3
      call yy_staggered_exchange(table,samples,value)
      if(repeat==1) saved=value
      if(any(saved/=value)) bad=bad+1
    end do
    err=0.d0;sumerr=0.d0
    do n=1,nt
      if(.not.table%covered(n)) then
        if(any(value(:,n)/=0.d0)) bad=bad+1
        cycle
      end if
      call reference_frame(points(2,n),points(3,n),patch,frame)
      xyz=points(1,n)*frame(:,1);call analytic(caseid,xyz,exact)
      err=max(err,maxval(abs(value(:,n)-exact)));sumerr=sumerr+sum((value(:,n)-exact)**2)
      write(u,'(3i9,3es26.17)') caseid,patch,indices(n),value(:,n)
    end do
    covered=count(table%covered)
    call MPI_ALLREDUCE(covered,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(err,peak,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(sumerr,globalerr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(caseid==0.and.peak>1.d-12) bad=bad+1
    if(rank==0) write(*,'(a,5i8,2es26.17)') 'RESULT ',scale,family,mode,caseid,total,peak,sqrt(globalerr/(3*total))
  end do
  close(u)
  ! Rebuild and exercise no-local-receiver ranks with nonempty donor traffic.
  call yy_staggered_clear(table)
  if(rank==0) then
    call yy_staggered_build(table,MPI_COMM_WORLD,patch,family,bn,nn,radial,angular,directory,points(:,:1))
  else
    call yy_staggered_build(table,MPI_COMM_WORLD,patch,family,bn,nn,radial,angular,directory,points(:,:0))
  end if
  deallocate(samples,value);allocate(samples(size(table%request,2)),value(3,size(table%covered)))
  do n=1,size(samples)
    samples(n)=dot_product([0.2d0,-0.3d0,0.5d0],table%donor_basis(:,n))
  end do
  call yy_staggered_exchange(table,samples,value)
  if(rank==0) then
    if(maxval(abs(value(:,1)-[0.2d0,-0.3d0,0.5d0]))>1.d-12) bad=bad+1
  end if
  call yy_staggered_clear(table)
  call MPI_ALLREDUCE(geom,geomall,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(geomall>1.d-13) bad=bad+1
  call MPI_ALLREDUCE(bad,badall,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(rank==0) write(*,'(a,i8,es26.17)') 'CHECK ',badall,geomall
  if(badall/=0) call MPI_ABORT(MPI_COMM_WORLD,22,ierr)
  call MPI_FINALIZE(ierr)
contains
  pure subroutine reference_frame(theta,phi,p,b)
    double precision,intent(in) :: theta,phi
    integer,intent(in) :: p
    double precision,intent(out) :: b(3,3)
    double precision :: rotation(3,3),local(3,3)
    local(:,1)=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
    local(:,2)=[cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)]
    local(:,3)=[-sin(phi),cos(phi),0.d0]
    rotation=0.d0;rotation(1,1)=-1.d0;rotation(2,3)=1.d0;rotation(3,2)=1.d0
    b=local
    if(p==1) b=matmul(rotation,local)
  end subroutine
  pure subroutine analytic(which,x,b)
    integer,intent(in) :: which
    double precision,intent(in) :: x(3)
    double precision,intent(out) :: b(3)
    b=[0.2d0,-0.3d0,0.5d0]
    if(which==1) b=[sin(x(2)),sin(x(3)),sin(x(1))]
  end subroutine
end program
