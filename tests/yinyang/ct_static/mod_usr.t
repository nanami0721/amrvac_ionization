! Actual experimental Yin-Yang UCT-HLL evolution on a static MHD background.
module mod_usr
  use mod_mhd
  use mod_physics, only: phys_to_conserved,phys_get_pthermal
  use mod_yinyang, only: yy_basis,yy_position
  use mod_functions_bfield, only: get_divb
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(state),allocatable :: initial_state(:)
  integer :: stage_calls=0,last_iteration=-1,history_unit,test_world_rank
  logical :: spherical_control=.false.
  double precision :: largest_emf=0.d0
contains
  subroutine usr_init()
    integer :: n
    namelist /static_test/ spherical_control
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,static_test,end=10)
10    close(unitpar)
    end do
    call MPI_COMM_RANK(MPI_COMM_WORLD,test_world_rank,ierrmpi)
    if(spherical_control)then
      call set_coordinate_system('spherical')
    else
      call set_coordinate_system('yinyang')
    end if
    call mhd_activate()
    usr_init_one_grid=>initial
    usr_init_vector_potential=>potential
    usr_ct_stage=>stage_diagnostic
    usr_process_global=>diagnostic
    usr_process_adv_global=>advanced
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
      call yy_basis(x(ix^D,2),x(ix^D,3),yy_patch,frame)
      call yy_position(x(ix^D,1),x(ix^D,2),x(ix^D,3),yy_patch,pos)
      ! A = B cross r / 2, B=(0.2,-0.1,0.3).
      cart=0.5d0*[-0.1d0*pos(3)-0.3d0*pos(2), &
        0.3d0*pos(1)-0.2d0*pos(3),0.2d0*pos(2)+0.1d0*pos(1)]
      if(iprob==0) cart=0.d0
      A(ix^D)=dot_product(cart,frame(:,idir))
    {end do\}
  end subroutine
  subroutine stage_diagnostic(qt,qdt,sCT,s,fluxes)
    double precision,intent(in) :: qt,qdt
    type(state),intent(in) :: sCT(max_blocks),s(max_blocks)
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    integer :: n,g
    stage_calls=stage_calls+1
    do n=1,igridstail
      g=igrids(n)
      if(.not.all(ieee_is_finite(fluxes(g)%edge))) call mpistop('nonfinite static-test EMF')
      largest_emf=max(largest_emf,maxval(abs(fluxes(g)%edge)))
    end do
  end subroutine
  subroutine advanced(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    call diagnostic(iit+1,qt+dt)
  end subroutine
  subroutine diagnostic(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    integer :: n,g,i,j,k,c,bc(3),gi(3),bn(3),nn(3),lo(3),hi(3),fl(3),u,callsmin,callsmax
    double precision :: divb(ixG^T),pth(ixG^T),local(11),maximum(11),mins(2),minimum(2),tot(3),total(3),error,rms,analytic,frame(3,3),expected(3)
    logical :: seam,radial,active
    character(len=80) :: name
    if(iit==last_iteration) return
    last_iteration=iit
    if(.not.allocated(initial_state)) then
      if(iit/=0)call mpistop('static test missing initial record')
      allocate(initial_state(max_blocks))
      do n=1,igridstail
        g=igrids(n);initial_state(g)%w=ps(g)%w;initial_state(g)%ws=ps(g)%ws
      end do
      if(test_world_rank==0)open(newunit=history_unit,file='static-history.txt',status='replace')
      call sleep(1)
    end if
    local=0.d0;mins=huge(1.d0);tot=0.d0
    lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3]
    bn=[block_nx1,block_nx2,block_nx3];nn=[domain_nx1,domain_nx2,domain_nx3]
    do n=1,igridstail
      g=igrids(n);block=>ps(g)
      bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
      if(.not.all(ieee_is_finite(ps(g)%w)))call mpistop('nonfinite static-test state')
      call get_divb(ps(g)%w,ixG^LL,ixG^LL,divb)
      call phys_get_pthermal(ps(g)%w,ps(g)%x,ixG^LL,ixM^LL,pth)
      do k=ixGlo3,ixGhi3
        do j=ixGlo2,ixGhi2
          do i=ixGlo1,ixGhi1
            gi=bc*bn+[i,j,k]-lo
            active=all([i,j,k]>=lo).and.all([i,j,k]<=hi)
            radial=gi(1)<0.or.gi(1)>=nn(1)
            if(radial)local(4)=max(local(4),abs(divb(i,j,k)))
            if(.not.active)then
              local(3)=max(local(3),abs(divb(i,j,k)))
              cycle
            end if
            seam=any(gi(2:3)<2).or.any(gi(2:3)>=nn(2:3)-2)
            local(1)=max(local(1),abs(divb(i,j,k)))
            if(seam)local(2)=max(local(2),abs(divb(i,j,k)))
            error=sqrt(sum((ps(g)%w(i,j,k,mag(:))-initial_state(g)%w(i,j,k,mag(:)))**2))
            local(5)=max(local(5),error)
            local(6)=max(local(6),sqrt(sum(ps(g)%w(i,j,k,mom(:))**2))/ps(g)%w(i,j,k,rho_))
            if(seam)local(7)=max(local(7),error)
            if(.not.seam.and.gi(1)>=2.and.gi(1)<nn(1)-2)local(8)=max(local(8),error)
            call yy_basis(ps(g)%x(i,j,k,2),ps(g)%x(i,j,k,3),yy_patch,frame)
            expected=[0.2d0,-0.1d0,0.3d0];if(iprob==0)expected=0.d0
            analytic=sqrt(sum((matmul(frame,ps(g)%w(i,j,k,mag(:)))-expected)**2))
            local(11)=max(local(11),analytic)
            tot=tot+[error**2*ps(g)%dvolume(i,j,k),ps(g)%dvolume(i,j,k),analytic**2*ps(g)%dvolume(i,j,k)]
            mins=min(mins,[ps(g)%w(i,j,k,rho_),pth(i,j,k)])
          end do
        end do
      end do
      do c=1,3
        fl=lo;fl(c)=fl(c)-1
        local(10)=max(local(10),maxval(abs(ps(g)%ws(fl(1):hi(1),fl(2):hi(2),fl(3):hi(3),c) &
          -initial_state(g)%ws(fl(1):hi(1),fl(2):hi(2),fl(3):hi(3),c))))
      end do
    end do
    local(9)=largest_emf
    call MPI_ALLREDUCE(local,maximum,11,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(mins,minimum,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(tot,total,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(stage_calls,callsmin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(stage_calls,callsmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(callsmin/=3*iit.or.callsmax/=3*iit)call mpistop('incorrect static-test RK stage count')
    rms=sqrt(total(1)/total(2))
    if(test_world_rank==0)then
      write(history_unit,'(i6,16es25.16)')iit,qt,maximum,rms,sqrt(total(3)/total(2)),minimum
      flush(history_unit)
      if(mod(iit,10)==0.or.iit==it_max)write(*,'(a,i6,16es25.16)')'CT_STATIC ',iit,qt,maximum,rms,sqrt(total(3)/total(2)),minimum
    end if
    if(maximum(1)>1.d-10.or.minval(minimum)<=0.d0)call mpistop('static-test divergence or positivity failure')
    if(iprob==0.and.maximum(5)>1.d-13)call mpistop('zero-field magnetic drift')
    if(iit==it_max)then
      write(name,'(a,i0,a)')'static-rank-',test_world_rank,'.txt'
      open(newunit=u,file=trim(name),status='replace')
      do n=1,igridstail
        g=igrids(n);bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
        do k=ixMlo3,ixMhi3
          do j=ixMlo2,ixMhi2
            do i=ixMlo1,ixMhi1
              gi=bc*bn+[i,j,k]-lo
              write(u,'(4i8,8es25.16)')yy_patch,gi,ps(g)%w(i,j,k,:)
            end do
          end do
        end do
      end do
      close(u)
    end if
  end subroutine
end module
