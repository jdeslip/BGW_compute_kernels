!==============================================================================
! kernel_BSE.x - Kernel for the BSE interpolation from BerkeleyGW
!                Kernel originally by Felipe H. da Jornada
!                (jornada@berkeley.edu)
!
! Run ./kernel_BSE.x for the documentation.
!
!==============================================================================

program kernel_BSE

  use mpi
  use common_m
  use intkernel_m
  use timing_m

  implicit none

  type(crystal) :: crys
  type(grid) :: kg
  type (xctinfo) :: xct
  complex(DPC), allocatable :: dcc(:,:,:,:,:), dvv(:,:,:,:,:)
  complex(DPC), allocatable :: hmtrx(:,:)
  complex(DPC) :: tmp(2,2)
  real(DP), allocatable :: intp_coefs(:,:)
  integer, allocatable :: fi2co_wfn(:,:)

  integer :: kgrid_co(3), kgrid_fi(3), ik, is, iv, ic, ikcvs, nkpt_own
  integer :: nrows, ncols
  integer :: i1, i2, i3, ii

  character*16, allocatable :: routnam(:)
  character*16 :: arg
  integer, allocatable :: routsrt(:)
  real(DP) :: tsec(2)
  integer :: ncount, narg, iargc
  integer :: nthreads, tid, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  integer :: nnodes, nranks, nranks_per_node, mpierr, inode

  inode = 0
  call MPI_Init(mpierr)
  if(mpierr .ne. MPI_SUCCESS) then
    write(0,'(a)') 'ERROR: MPI initialization failed!'
    stop 999
  endif
  call MPI_Comm_size(MPI_COMM_WORLD, nranks_per_node, mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, inode, mpierr)

  if (inode==0) then
    write(6,*)
    write(6,*) 'BerkeleyGW - Kernel for the interpolation of the BSE matrix'
    write(6,*) '==========================================================='
    write(6,*)
  endif
  narg = iargc()
  if (narg/=4 .and. narg/=7) then
    if (inode==0) then
      write(6,*) 'Usage:'
      write(6,*) ' (1) kernel_BSE.x  nnodes  nvb_co ncb_co nkx_co  nvb_fi ncb_fi nkx_fi; or'
      write(6,*) ' (2) kernel_BSE.x  nnodes  nvb ncb nkx'
      write(6,*)
      write(6,*) 'Short examples (<1s):'
      write(6,*) ' (S1) ./kernel_BSE.x  1  2 2 2  1 1 4'
      write(6,*) ' (S2) ./kernel_BSE.x  1  10 25 2'
      write(6,*) ' (S3) ./kernel_BSE.x  1  30 60 1'
      write(6,*)
      write(6,*) 'Long examples (<30s):'
      write(6,*) ' (L1) ./kernel_BSE.x  1  6 15 2  2 8 4'  ! last param is sqrt(nkpt_fi)
      write(6,*) ' (L1.1) ./kernel_BSE.x  1  6 15 4  2 8 10'  ! last param is sqrt(nkpt_fi) 3rd param is sqrt(nkpt_co)
      write(6,*) ' (L2) ./kernel_BSE.x  1  2 2 3  1 1 7'
      write(6,*) ' (L3) ./kernel_BSE.x  1  8 12 5'
      write(6,*)
      write(6,*) 'These examples were prioritized, so (S1)/(L1) are more relevant than (S3)/(L3).'
      write(6,*)
      call flush(6)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    stop
  endif

  ! FHJ - Kernel initialization
  ! ---------------------------------------------------------------------------

  xct%skipinterp = .true.
  xct%npts_intp_kernel = 1
  call getarg(1, arg)
  read(arg,*) nnodes
  call getarg(2, arg)
  read(arg,*) xct%nvb_co
  call getarg(3, arg)
  read(arg,*) xct%ncb_co
  call getarg(4, arg)
  read(arg,*) kgrid_co(1)
  kgrid_co(2) = kgrid_co(1)
  kgrid_co(3) = kgrid_co(1)
  xct%nvb_fi = xct%nvb_co
  xct%ncb_fi = xct%ncb_co
  kgrid_fi(:) = kgrid_co(:)
  if (narg==7) then
    xct%skipinterp = .false.
    xct%npts_intp_kernel = 3
    call getarg(5, arg)
    read(arg,*) xct%nvb_fi
    call getarg(6, arg)
    read(arg,*) xct%ncb_fi
    call getarg(7, arg)
    read(arg,*) kgrid_fi(1)
    kgrid_fi(2) = kgrid_fi(1)
    kgrid_fi(3) = kgrid_fi(1)
  endif

  xct%nkpt_co = product(kgrid_co)
  xct%nkpt_fi = product(kgrid_fi)
  xct%nktotal = xct%nkpt_fi
  xct%nspin = 1
  xct%n1b_co = xct%nvb_co
  xct%n2b_co = xct%ncb_co
  xct%idimensions = 3
  nranks = nranks_per_node*nnodes

  nkpt_own = (xct%nkpt_fi+nranks-1)/nranks
  nrows = xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin
  ncols = nkpt_own*xct%ncb_fi*xct%nvb_fi*xct%nspin

  if (inode==0) then
    write(6,*) 'Input parameters'
    write(6,*) '----------------'
    write(6,'(1x,a," = ",l1)') 'skip_interpolation', xct%skipinterp
    write(6,666) 'nvb_co', xct%nvb_co
    write(6,666) 'ncb_co', xct%ncb_co
    write(6,666) 'kgrid_co(1)', kgrid_co(1)
    write(6,666) 'nkpt_co', xct%nkpt_co
    write(6,666) 'nvb_fi', xct%nvb_fi
    write(6,666) 'ncb_fi', xct%ncb_fi
    write(6,666) 'kgrid_fi(1)', kgrid_fi(1)
    write(6,666) 'nkpt_fi', xct%nkpt_fi
    write(6,666) 'nkpt_own', nkpt_own
    write(6,666) 'nrows', nrows
    write(6,666) 'ncols', ncols
    write(6,666) 'nnodes', nnodes
    write(6,666) 'nranks_per_node', nranks_per_node
    write(6,666) 'nranks', nranks
  endif
!$omp parallel private(tid) shared(nthreads)
  tid = OMP_GET_THREAD_NUM()
  if (tid==0) then
    nthreads = OMP_GET_NUM_THREADS()
    if (inode==0) then
      write(6,666) 'Number of threads', nthreads
    endif
  endif
!$omp end parallel
  666 format(1x,a,' = ',i0)
  if (inode==0) then
    write(6,*)
  endif

  call timacc(0,0)
  call timacc(1,1)

  if (inode==0) then
    write(6,*) 'Initializing default params'
  endif
  call timacc(2,1)

  crys%celvol = 1d0
  crys%bdot = 0d0
  crys%bdot(1,1) = 1d0
  crys%bdot(2,2) = 1d0
  crys%bdot(3,3) = 1d0
  xct%wcoul0 = 10d0

  allocate(hmtrx(nrows,ncols))
  if (inode==0) then
    write(6,*)"Size of hmtrx/nkpt_fi = ",(16.0*nrows*ncols/(xct%nkpt_fi))/1024.0," kbytes"
  endif

  hmtrx(:,:) = (0d0, 0d0)
  ! FHJ: Note that the "distributed" hmtrx matrix makes no sense!
  do ik=1, nkpt_own
    do is=1,xct%nspin
      do iv=1,xct%nvb_fi
        do ic=1,xct%ncb_fi
          ikcvs = bse_index(ik, ic, iv, is, xct)
          hmtrx(ikcvs,ikcvs) = 1d0 + ic*1d0 + iv*1d0 + (ik-xct%nkpt_fi*0.5d0)*1d-3
        enddo
      enddo
    enddo
  enddo

  allocate(kg%f(3,xct%nkpt_fi))
  allocate(fi2co_wfn(xct%npts_intp_kernel,xct%nkpt_fi))
  allocate(intp_coefs(xct%npts_intp_kernel, xct%nkpt_fi))
  ik=0
  do i1=0,kgrid_fi(1)-1
    do i2=0,kgrid_fi(2)-1
      do i3=0,kgrid_fi(3)-1
        ik=ik+1
        kg%f(1,ik)=(dble(i1))/dble(kgrid_fi(1))
        kg%f(2,ik)=(dble(i2))/dble(kgrid_fi(2))
        kg%f(3,ik)=(dble(i3))/dble(kgrid_fi(3))
        fi2co_wfn(1,ik) = int(ik*dble(xct%nkpt_co)/dble(xct%nkpt_fi)) + 1
        do ii=1, xct%npts_intp_kernel
          intp_coefs(ii,ik) = 1d0 / xct%npts_intp_kernel
          fi2co_wfn(ii,ik) = fi2co_wfn(1,ik) + ii-1
          fi2co_wfn(ii,ik) = mod(fi2co_wfn(ii,ik)-1, xct%nkpt_co) + 1
        enddo
      enddo
    enddo
  enddo

  allocate(dvv(xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
  allocate(dcc(xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
  do ii = 1, xct%npts_intp_kernel
    do ik = 1, xct%nkpt_fi
      do is = 1, xct%nspin
        do iv = 1, xct%nvb_fi
          dvv(iv, :, is, ik, ii) = ii*1d-1 + ik*1d-2 + is*1d-3 + iv*(1d-3,1d-2)
          dvv(iv, 1, is, ik, ii) = 1d0 - sum(abs(dvv(iv, :, is, ik, ii))**2)
        enddo
        do ic = 1, xct%ncb_fi
          dcc(ic, :, is, ik, ii) = ii*1d-1 + ik*1d-2 + is*1d-3 + ic*(1d-2,1d-3)
          dcc(ic, 1, is, ik, ii) = 1d0 - sum(abs(dcc(ic, :, is, ik, ii))**2)
        enddo
      enddo
    enddo
  enddo
  call timacc(2,2)


  ! FHJ - This is where we do the actual computation
  ! ---------------------------------------------------------------------------

  if (inode==0) then
    write(6,*) 'Calling intkernel'
  endif
  call timacc(5,1)
  call intkernel(crys,kg,xct,hmtrx,dcc,dvv,fi2co_wfn,intp_coefs,nkpt_own,inode)
  call timacc(5,2)
  if (inode==0) then
    write(6,*) 'Done!'
  endif


  ! FHJ - Summarize results
  ! ---------------------------------------------------------------------------

  call timacc(6,1)
  if (inode==0) then
    write(6,*)
    write(6,*) 'Summary of resulting matrix:'
    write(6,*) '----------------------------'
    write(6,*) 'Sum(H) = ', sum(hmtrx)
    write(6,*) 'Sum(|real(H)|) = ', sum(abs(dble(hmtrx)))
    write(6,*) 'Sum(|imag(H)|) = ', sum(abs(aimag(hmtrx)))
    write(6,*) 'Norm2(real(H)+imag(H)) = ', norm2(dble(hmtrx)+aimag(hmtrx))
    write(6,*) 'H(1:2,1:2) = '
    tmp = hmtrx(1:2, 1:2)
    write(6,*) tmp
    write(6,*)
  endif
  call timacc(6,2)


  ! FHJ - Timing
  ! ---------------------------------------------------------------------------

  allocate(routnam(60))
  routnam(1)='TOTAL:'
  routnam(2)='SETUP:'
  routnam(5)='INTKERNEL:'
  routnam(6)='SUMMARIZE:'
  routnam(51)='IK Setup:'
  routnam(56)='IK Cells:'
  routnam(57)='IK Interp:'
  routnam(58)='IK Sum:'
  allocate(routsrt(8))
  routsrt=(/ 2, 5, 6, 51, 56, 57, 58, 1 /)

  call timacc(1,2)
  if (inode==0) then
    write(6,*)
    write(6,9000) 'CPU (s)','WALL (s)','#'
    write(6,*)
    do ii=1,ubound(routsrt, 1)
      call timacc(routsrt(ii),3,tsec,ncount)
      if (ii>1) then
        if (abs(routsrt(ii)-routsrt(ii-1))>10) write(6,*)
      endif
      write(6,9001) routnam(routsrt(ii)),tsec(1),tsec(2),ncount
    enddo
    write(6,*)
  endif
  call flush(6)
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

9000 format(22x,a13,  3x,a13,  3x,a8)
9001 format(1x,a16,5x,f13.3,3x,f13.3,3x,i8)

end program kernel_BSE
