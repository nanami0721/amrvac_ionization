module mod_usr
  use mod_rho
  use mod_yinyang, only: yy_build_staggered_table, yy_basis, yy_position
  use mod_yinyang_staggered, only: yy_staggered_table, yy_staggered_offset, &
    yy_staggered_exchange, yy_staggered_clear
  implicit none
  logical :: checked=.false.
contains
  subroutine usr_init()
    call set_coordinate_system('yinyang')
    call rho_activate()
    usr_init_one_grid=>initial
    usr_process_global=>check_tables
  end subroutine
  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    w(ixO^S,rho_)=1.d0
  end subroutine
  subroutine check_tables(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    type(yy_staggered_table) :: table
    double precision,allocatable :: points(:,:),samples(:),values(:,:)
    double precision :: offset(3),coord(3),frame(3,3),xyz(3),error,maximum
    integer :: family,n,g,i,j,k,c,t,side,axis,idx(3),bn(3),nb(3),pig(3),lo(3),hi(3),countall,repeat
    if(checked) return
    checked=.true.;bn=[block_nx1,block_nx2,block_nx3];nb=[ng1(1),ng2(1),ng3(1)]
    lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3]
    allocate(points(3,12*igridstail))
    do family=1,2
      t=0
      do n=1,igridstail
        g=igrids(n);pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
        do axis=2,3
          do side=1,2
            if(side==1.and.pig(axis)/=1) cycle
            if(side==2.and.pig(axis)/=nb(axis)) cycle
            do c=1,3
              call yy_staggered_offset(family,c,offset)
              idx=lo+bn/2
              idx(axis)=lo(axis)-1
              if(side==2) idx(axis)=hi(axis)+1
              i=idx(1);j=idx(2);k=idx(3)
              coord=ps(g)%x(i,j,k,:)+(offset-0.5d0)*ps(g)%dx(i,j,k,:)
              t=t+1;points(:,t)=coord
            end do
          end do
        end do
      end do
      call yy_build_staggered_table(table,ps,family,points(:,:t))
      if(any(.not.table%covered)) call mpistop('actual AMRVAC staggered receiver lacks donors')
      allocate(samples(size(table%request,2)),values(3,t));error=0.d0
      do n=1,size(samples)
        g=table%request(1,n);idx=table%request(2:4,n)+lo;c=table%request(5,n)
        i=idx(1);j=idx(2);k=idx(3)
        call yy_staggered_offset(family,c,offset)
        coord=ps(g)%x(i,j,k,:)+(offset-0.5d0)*ps(g)%dx(i,j,k,:)
        call yy_basis(coord(2),coord(3),yy_patch,frame)
        call yy_position(coord(1),coord(2),coord(3),yy_patch,xyz)
        error=max(error,maxval(abs(xyz-table%donor_position(:,n))), &
          maxval(abs(frame(:,c)-table%donor_basis(:,n))))
        samples(n)=dot_product([0.2d0,-0.3d0,0.5d0],frame(:,c))
      end do
      do repeat=1,3
        call yy_staggered_exchange(table,samples,values)
        do n=1,t
          error=max(error,maxval(abs(values(:,n)-[0.2d0,-0.3d0,0.5d0])))
        end do
      end do
      call MPI_ALLREDUCE(error,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
      call MPI_ALLREDUCE(t,countall,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
      if(yy_world_rank==0) write(*,'(a,2i8,es26.17)') 'STAGGERED_MESH ',family,countall,maximum
      if(maximum>1.d-12) call mpistop('actual AMRVAC staggered table disagrees with mesh')
      deallocate(samples,values);call yy_staggered_clear(table)
    end do
  end subroutine
end module
