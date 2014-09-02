module intkernel_m

  use common_m
  use blas_m
  use cells_m
  use timing_m

  implicit none

  private
  
  public :: intkernel

contains

subroutine intkernel(crys,kg,xct,hmtrx,dcc,dvv,fi2co_wfn,intp_coefs)
  type (crystal), intent(in) :: crys
  type (grid), intent(in) :: kg
  type (xctinfo), intent(inout) :: xct
  complex(DPC), intent(inout) :: hmtrx(:,:) !< (xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin,peinf%nblocks*peinf%nblockd)
  !> (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel) Band expansion coefficients
  complex(DPC), intent(in) :: dcc(:,:,:,:,:), dvv(:,:,:,:,:)
  !> (xct%npts_intp_kernel, xct%nkpt_fi) For a given k-point in the fine grid
  !! returns the closest k-point in the coarse WFN k-grid.
  integer, intent(in) :: fi2co_wfn(:,:)
  !> (xct%npts_intp_kernel, xct%nkpt_fi) Delaunay/greedy interpolation coefficients
  real(DP), intent(in) :: intp_coefs(:,:)

  real(DP) :: vcoul, oneoverq

! FHJ: these arrays depend only on the distance w.r.t qq=0 (inside the BZ)
  real(DP), allocatable :: dist_array(:) !length of dq for a given index
  real(DP), allocatable :: vcoul_array(:), oneoverq_array(:)

  !> FHJ: cells structure that keeps all fine (k-kp) transitions
  type(cells_t) :: cells_fi
  ! FHJ: these variables store the cached value of |q|^2, v(q), epsinv, for 
  ! each dq in the same order as it appears in the cell structure.
  real(DP) :: dq(3), dq_bz(3), abs_q2
  real(DP), allocatable :: dqs(:,:), dqs_bz(:,:), abs_qs2(:), eps_intp(:)
  integer, allocatable :: vq_map(:)
  integer :: ik_cells, jkp_offset
  integer, allocatable :: ikb(:)
  
  integer :: icb,ivb,imatrix, &
    ik,ic,iv,ikp,icp,ivp,ikt,ikcvs,ikcvsd,jj,jc,jv,js, &
    jcp,jvp,jsp,jk,jkp,dimbse,nmatrices, &
    iold,icout,ivout,ncdim,nvdim,inew,ibt
  !> (xct%nkpt_co) maps a k-point in the WFN_co to the bse*mat files.
  integer, allocatable :: wfn2bse(:)
  !> (xct%npts_intp_kernel, xct%nkpt_fi) For a given k-point in the fine grid,
  !! returns the closest k-point in the (coarse) bse*mat files.
  integer, allocatable :: fi2co_bse(:,:)
  integer :: ivertp, ivert
  integer, allocatable :: jkp2offset(:)
  real(DP) :: bsemat_fac,fac_d,fac_x,w_eff,eps,factor

!------------------------------
! kernel and transformation matrices
  complex(DPC), allocatable :: &
    bsedmatrix(:,:,:,:,:,:,:),bsedmt(:,:,:,:,:,:,:), &
    bsedmatrix_loc(:,:,:,:,:,:,:), &
    dcckp(:,:,:), dvvkp(:,:,:)

  allocate(ikb(xct%nkpt_fi))
  ikb = (/ (ik, ik=1,xct%nkpt_fi) /)
  ibt = xct%nkpt_fi
  allocate(wfn2bse(xct%nkpt_co))
  wfn2bse = (/ (ik, ik=1,xct%nkpt_co) /)
  allocate(fi2co_bse(xct%npts_intp_kernel, xct%nkpt_fi))

  factor = -8.d0*PI_D/(crys%celvol*xct%nktotal)
  fac_d = 1d0
  fac_x = -1d0
  
!-----------------------------------------------------------------------------!
!              EPSILON INTERPOLATION AND V(Q) PRE-CALCULATION
! FHJ: The interpolation of epsinv scales as O(nkpt) for regular grids.
! We pre-calculate all dqs = k(:,:) - k(:,1), move this quantity to the [0, 1) 
! range, and put them into a cell structure. For each dq, we interpolate epsinv.
! In a separate memory copy, we move all dqs to the 1st BZ (dqs_bz), and
! calculate v(dq_bz) and related quantities.
! When we reach the actual kernel interpolation, we calculate dq=k-kp, move it
! to the same [0, 1) range, and get the mapping to the dqs, which is done in
! O(1) using the cell structure.
! This algorithm should also works for coarse-grained grids (a fine grid with 
! holes), as long as the dqs are compute with all possible fine transitions.

  allocate(eps_intp(xct%nktotal))
  allocate(dist_array(xct%nktotal))
  allocate(vcoul_array(xct%nktotal))
  allocate(oneoverq_array(xct%nktotal))
  allocate(dqs(3, xct%nktotal))
  allocate(dqs_bz(3, xct%nktotal))
  allocate(abs_qs2(xct%nktotal))
  allocate(vq_map(xct%nktotal))

  ! FHJ: This block is equivalent to dqs(:,ik) = UNIT_RANGE(kg%f(:,ik)-kg%f(:,1))
  ! if nkpt_fi==nktotal, but only up to TOL_SMALL. In order for the testsuite not to
  ! fail, we explicitly break this into standard and patched_sampling cases.
  do ik = 1, xct%nkpt_fi 
    dqs(:,ik) = kg%f(:,ik) - kg%f(:,1)
    dqs(:,ik) = dqs(:,ik) - floor(dqs(:,ik))
  enddo

  ! FHJ: TODO: parallelize this loop (trivial)
  call timacc(54,1)
  call cells_init(cells_fi, dqs, periodic=.true., quiet=.true.)
  do ik = 1, xct%nktotal
    dq(:) = dqs(:, ik)
    
    ! FHJ: Move pt to 1st BZ and save it
    call point_to_bz(crys, dq, dq_bz, abs_q2)
    dqs_bz(:, ik) = dq_bz(:)
    abs_qs2(ik) = abs_q2

  ! Interpolate epsinv - Just the head
    eps_intp(ik) = exp(-abs_q2*1d-2)
  enddo !ik
  call timacc(54,2)
! FHJ: Done with epsinv interpolation
!-----------------------------------

! FHJ: v(q) pre-calculation
!-------------------------------
! We should all contribute to calculating vcoul0 since it might involve minibz average

  call timacc(55,1)
  inew = 0

  do ik = 1, xct%nktotal
    dq(:) = dqs(:, ik)
    dq_bz(:) = dqs_bz(:, ik)
    abs_q2 = abs_qs2(ik)

    ! FHJ: Take a look if v(q) for this |q| has already been calculated
    ! We assume that v(q) = v(|q|), which is fine since we only include q`s
    ! in the periodic directions, so truncation doesn`t affect.
    iold=0
    do jj=inew,1,-1
      if (dabs(dist_array(jj)-abs_q2) < TOL_Small) then
        iold=jj
        exit
      endif
    enddo

    if (iold==0) then
      inew=inew+1
      dist_array(inew)=abs_q2
      oneoverq = 1d0/sqrt(abs_q2+TOL_SMALL)
      vcoul_array(inew) = 1d0/(abs_q2+TOL_SMALL)
      oneoverq_array(inew) = oneoverq/(8.0*Pi_D)
      vq_map(ik) = inew
    else
      vq_map(ik) = iold
    endif !iold == 0
  enddo !ik

! FHJ: Done with v(q) pre-calculation
!-----------------------------------
  write(6,'(1x,a,i6,a)') 'Finished calculating Vcoul with ',inew,' unique points'

  call timacc(55,2)
  call timacc(51,1)

!--------------------------------
! Allocate data

  if (.not. xct%skipinterp) then
    write(6,'(1x,a)') 'Performing Kernel Interpolation'
  else
    write(6,'(1x,a)') 'Building Interaction Kernel'
  endif

!  if (xct%ipar .eq. 1) then
    allocate(bsedmt(xct%nvb_fi,xct%ncb_fi,xct%nvb_fi,xct%ncb_fi,xct%nspin,xct%nspin,1))
!  else if (xct%ipar .eq. 2) then
!    allocate(bsedmt(xct%nvb_fi,xct%ncb_fi,xct%nvb_fi,1,xct%nspin,xct%nspin,1))
!  else
!    allocate(bsedmt(xct%nvb_fi,xct%ncb_fi,1,1,xct%nspin,xct%nspin,1))
!  endif
  allocate(dcckp(xct%n2b_co,xct%ncb_fi,xct%nspin))
  allocate(dvvkp(xct%n1b_co,xct%nvb_fi,xct%nspin))
 
  dimbse=xct%nkpt_co*(xct%n2b_co*xct%n1b_co*xct%nspin)**2

  call timacc(51,2)
  call timacc(52,1)

  do ik = 1, xct%nkpt_fi
    fi2co_bse(:,ik) = wfn2bse(fi2co_wfn(:,ik))
  enddo


  call timacc(52,2)

!------------- Read in coarse matrices: head, wings, body, exchange --------------------
  nmatrices = 4
  write(6,'(1x,a,i0,a)') 'Interpolating BSE kernel with ', ibt*nmatrices, ' blocks'

  ! FHJ: For each PE, figure our which coarse k-points jkp it needs in order to
  ! interpolate the BSE matrix on the fine k-points it owns.
  allocate(jkp2offset(xct%nkpt_co))
  jkp2offset = -1 ! offset=-1 means we don`t need this coarse k-point
  jkp_offset = 0
  do jkp=1,xct%nkpt_co
    jkp2offset(jkp) = jkp_offset
    jkp_offset = jkp_offset + xct%nkpt_co 
  enddo

  allocate(bsedmatrix_loc(xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nspin,xct%nspin,jkp_offset))
  nmatrices = 4
  do imatrix = 1, nmatrices     
    allocate(bsedmatrix(xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nspin,xct%nspin,xct%nkpt_co))
    do jkp=1,xct%nkpt_co

          do jj=1, xct%n2b_co * xct%n1b_co
            ! FHJ: we would normally be reading a file here..
            jvp = mod(jj-1, xct%n1b_co) + 1
            jcp = (jj-1)/xct%n1b_co + 1
            ikp = jkp
            do jk = 1, xct%nkpt_co
              do jc = 1, xct%n2b_co
                do jv = 1, xct%n1b_co
                  do js = 1, xct%nspin
                    do jsp = 1, xct%nspin
                      bsedmatrix(jv,jc,jvp,jcp,js,jsp,jk) = cmplx( (js+jsp) + 2d0*abs(jv-jvp), 2d0*(jc-jcp) + 16d0*(jkp-jk))
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
      
      if (jkp2offset(jkp)/=-1) then
        bsedmatrix_loc(:,:,:,:,:,:,jkp2offset(jkp)+1:jkp2offset(jkp)+xct%nkpt_co) = &
            bsedmatrix(:,:,:,:,:,:,1:xct%nkpt_co)
      endif
    enddo
    deallocate(bsedmatrix)

!==============================================================================
! FHJ: Beginning of the big loop over k-points
!      For each fine k-points ik/ikp and expansion vertices ivert/iterp, find
!      the corresponding coarse-grid points jk/jkp and interpolate the kernel.
!==============================================================================
    do ikt = 1, ibt
      ikp = ikb(ikt)
      do ivertp = 1, xct%npts_intp_kernel
        jkp = fi2co_bse(ivertp, ikp)
        jkp_offset = jkp2offset(jkp)

          do js=1,xct%nspin
            dcckp(:,:,js) = transpose(dcc(:,:,js,ikp,ivertp))
            dvvkp(:,:,js) = transpose(dvv(:,:,js,ikp,ivertp))
          enddo

            ! FHJ: please don`t fix the indentation yet
            do ik = 1, xct%nkpt_fi
              do ivert = 1, xct%npts_intp_kernel
                jk = fi2co_bse(ivert, ik)
              
                bsedmt = (0d0,0d0)

                if (ivert==1) then
                ! FHJ: No need to update this block for different jk`s
                dq(:) = kg%f(:,ik) - kg%f(:,ikp)

                ! FHJ: and now we can reuse all the information from the cell
                ! structure to get the mapping to the BZ, v(q) and the interpolated
                ! epsilon in O(1) operations.
                call timacc(56,1)
                call cells_find_exactly(cells_fi, dq, ik_cells)
                if (ik_cells==0) then
                  write(0,'(a)') 'Found a point that was not calculated before:'
                  write(0,'(3(f13.9,1x))') dq
                  write(0,'(a)') 'Are you using a non-uniform grid?'
                  stop
                endif
                abs_q2 = abs_qs2(ik_cells)
                eps = eps_intp(ik_cells)
                iold = vq_map(ik_cells)

!------------------- Calculate Coulomb Interaction-------------------------------------

                ! The 1/(8*pi) is already included
                vcoul = vcoul_array(iold)
                oneoverq = oneoverq_array(iold)
                call timacc(56,2)
                
                w_eff = vcoul * eps
                if (abs_q2<TOL_Zero) then
                  w_eff = xct%wcoul0/(8.0*PI_D)
                endif

!---------------- Start interpolation for the present pair (k,kp)----------------------------

!-----------------------------
! Head (spin diagonal)

!                if (xct%ipar .eq. 1) then
                  ncdim=xct%ncb_fi
                  nvdim=xct%nvb_fi
                  icout=1
                  ivout=1
!                else if (xct%ipar .eq. 2) then
!                  ncdim=1
!                  nvdim=xct%nvb_fi
!                  icout=peinf%icb(peinf%inode+1,ikt)
!                  ivout=1
!                else if (xct%ipar .eq. 3) then
!                  ncdim=1
!                  nvdim=1
!                  icout=peinf%icb(peinf%inode+1,ikt)
!                  ivout=peinf%ivb(peinf%inode+1,ikt)
!                endif
                endif !ivert==1

                call timacc(57,1)

                !write(6,*) bsedmatrix_loc(1,1,1,1,1,1,jkp_offset+1)
                if (.not. xct%skipinterp) then
                  call interpolate(xct, &
                    bsedmatrix_loc(:,:,:,:,:,:,jkp_offset+jk), bsedmt(:,:,:,:,:,:,1), &
                    dcc(:,:,:,ik,ivert), dcckp, dvv(:,:,:,ik,ivert), dvvkp)
                else
! XXX Thread?
                    do jsp=1,xct%nspin
                      do js=1,xct%nspin
                        do jcp=1,ncdim
                          do jvp=1,nvdim
                            do jc=1,xct%ncb_fi
                              do jv=1,xct%nvb_fi
                                bsedmt(jv,jc,jvp,jcp,js,jsp,1) = bsedmatrix_loc( &
                                  jv, jc, jvp+ivout-1, jcp+icout-1, &
                                  js, jsp, jkp_offset+jk)
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                endif

                call timacc(57,2)

!------------------------------
! Add interaction kernel to the Hamiltonian

                if (imatrix .eq. 1) then
                    ! change back
                    bsemat_fac = fac_d * w_eff
                else if (imatrix .eq. 2) then
                    bsemat_fac = fac_d * oneoverq
                else if (imatrix .eq. 3) then
                  bsemat_fac = fac_d
                else if (imatrix .eq. 4) then
                  bsemat_fac = -fac_x
                  if (xct%nspin .eq. 1) bsemat_fac = bsemat_fac * 2D0
                endif
                bsemat_fac = bsemat_fac * intp_coefs(ivert, ik) * intp_coefs(ivertp, ikp)

                if(abs(bsemat_fac) < TOL_Zero) cycle
                call timacc(58,1)

! JRD: When comparing these loops with the loops in diag.f90 notice that:
! iv -> ivp
! ic -> icp
! ...
! 
! hmtrx array is accessed the same way since ikcvs here is based on ik,ic,iv and
! ikcvsd is based on ikp,icp,ivp. This is the opposite of diag.f90 which makes the 
! loops identical. Someone should fix the variable names to match at some point...

! XXX: Thread?
! XXX: WHY IS BSE_INDEX FASTER ON THE SPIN INDEX??

                do icp = 1, xct%ncb_fi
!                  if (xct%ipar .eq. 1 .or. icp .eq. 1) then
                    do ivp = 1, xct%nvb_fi
                      if (xct%ipar .lt. 3 .or. ivp .eq. 1) then
                        do jsp = 1, xct%nspin
!                          if (xct%ipar .eq. 1) then 
                            icb=icp
                            ivb=ivp
                            ikcvsd = bse_index(ikt, icb, ivb, jsp, xct)
!                          else if (xct%ipar .eq. 2) then
!                            !icb=peinf%icb(peinf%inode+1,ikt)
!                            icb=1
!                            ivb=ivp
!                            ikcvsd = bse_index(ikt, icb, ivb, jsp, xct, ncband = 1)
!                          else
!                            icb=1
!                            ivb=1
!                            ikcvsd = bse_index(ikt, icb, ivb, jsp, xct, ncband = 1, nvband = 1)
!                          endif

                          do ic = 1, xct%ncb_fi
                            do iv = 1, xct%nvb_fi
                              do js = 1,xct%nspin
                                ikcvs = bse_index(ik, ic, iv, js, xct)
                                hmtrx(ikcvs,ikcvsd) = hmtrx(ikcvs,ikcvsd) + &
                                  bsemat_fac * bsedmt(iv,ic,ivb,icb,js,jsp,1)
                              enddo              ! js
                            enddo              ! iv
                          enddo              ! ic
                        enddo              ! jsp
                      endif
                    enddo              ! ivp
!                  endif
                enddo              ! icp
                call timacc(58,2)
              enddo ! ivert / jk
            enddo ! ik
      enddo ! ivertp / jkp
    enddo ! ikp
  enddo ! imatrix 
!==============================================================================
! FHJ: End of the big loop over k-points
!==============================================================================

  call cells_free(cells_fi)

  deallocate(bsedmatrix_loc)
  deallocate(bsedmt)
  deallocate(dcckp)
  deallocate(dvvkp)
  deallocate(dist_array)
  deallocate(vcoul_array)
  deallocate(oneoverq_array)
  deallocate(dqs)
  deallocate(dqs_bz)
  deallocate(abs_qs2)
  deallocate(eps_intp)
  deallocate(vq_map)
  deallocate(fi2co_bse)
  deallocate(jkp2offset)

  return
end subroutine intkernel

!=======================================================================================

subroutine interpolate(xct,bse_co,bse_fi,dcck,dcckp,dvvk,dvvkp)
  type (xctinfo), intent(inout) :: xct
     ! xct should not be changed but ifort 12.0.0 -O3 compiles
     ! incorrectly if it is intent(in)
  complex(DPC), intent(in) :: bse_co(:,:,:,:,:,:) !< (xct%n1b_co,xct%n2b_co,xct%n1b_co,xct%n2b_co,xct%nspin,xct%nspin)
  complex(DPC), intent(out) :: bse_fi(:,:,:,:,:,:) !< (xct%nvb_fi,xct%ncb_fi,nvdim,ncdim,xct%nspin,xct%nspin)
  complex(DPC), intent(in) :: dcck(:,:,:)  !< (xct%ncb_fi,xct%n2b_co,xct%nspin)
  complex(DPC), intent(in) :: dcckp(:,:,:) !< (xct%n2b_co,xct%ncb_fi,xct%nspin)
  complex(DPC), intent(in) :: dvvk(:,:,:)  !< (xct%nvb_fi,xct%n1b_co,xct%nspin)
  complex(DPC), intent(in) :: dvvkp(:,:,:) !< (xct%n1b_co,xct%nvb_fi,xct%nspin)

  integer :: js,jsp,iv,ivp,jc,jcp,&
    js_dvvk, js_dvvkp, js_dcck, js_dcckp, bse_co_js, bse_co_jsp
  
  complex(DPC), allocatable :: mat_vcvc(:,:,:,:),mat_vfvc(:,:,:,:), &
    mat_vfvf(:,:,:,:),mat_cccc(:,:,:,:), &
    mat_cfcc(:,:,:,:),mat_cfcf(:,:,:,:)
!  complex(DPC), allocatable :: dummy(:,:,:,:), dummyp(:,:,:,:), dummy2(:,:,:,:), dummy3(:,:,:,:), dvvkn(:,:)
  
  bse_fi=0.0
  
  do js=1,xct%nspin
    do jsp=1,xct%nspin
      
      if (xct%nspin .eq. 1) then
        js_dcck=1
        js_dcckp=1
        js_dvvk=1
        js_dvvkp=1
        bse_co_js=1
        bse_co_jsp=1
      end if

! Fastest but worst on memory

        allocate(mat_vcvc(xct%n1b_co,xct%n1b_co,xct%n2b_co,xct%n2b_co))

        do jcp=1,xct%n2b_co
          do jc=1,xct%n2b_co
            mat_vcvc(:xct%n1b_co,1:xct%n1b_co,jc,jcp) = &
              bse_co(1:xct%n1b_co,jc,1:xct%n1b_co,jcp,bse_co_js,bse_co_jsp)
          enddo
        enddo

! Interpolate v,v`

        allocate(mat_vfvc(xct%nvb_fi,xct%n1b_co,xct%n2b_co,xct%n2b_co))

        do jc=1,xct%n2b_co
          do jcp=1,xct%n2b_co
            mat_vfvc(1:xct%nvb_fi,1:xct%n1b_co,jc,jcp) = &
              MATMUL((dvvk(1:xct%nvb_fi,1:xct%n1b_co,js_dvvk)), &
              (mat_vcvc(1:xct%n1b_co,1:xct%n1b_co,jc,jcp)))
          enddo
        enddo

        deallocate(mat_vcvc)
        allocate(mat_vfvf(xct%nvb_fi,xct%nvb_fi,xct%n2b_co,xct%n2b_co))

        do jc=1,xct%n2b_co
          do jcp=1,xct%n2b_co
            mat_vfvf(1:xct%nvb_fi,1:xct%nvb_fi,jc,jcp) = &
              MATMUL((mat_vfvc(1:xct%nvb_fi,1:xct%n1b_co,jc,jcp)), &
              conjg(dvvkp(1:xct%n1b_co,1:xct%nvb_fi,js_dvvkp)))
          enddo
        enddo

! Reorder from v,v` to c,c`

        deallocate(mat_vfvc)
        allocate(mat_cccc(xct%n2b_co,xct%n2b_co,xct%nvb_fi,xct%nvb_fi))

        do jcp=1,xct%n2b_co
          do jc=1,xct%n2b_co
            mat_cccc(jc,jcp,1:xct%nvb_fi,1:xct%nvb_fi) = &
              mat_vfvf(1:xct%nvb_fi,1:xct%nvb_fi,jc,jcp)
          enddo
        enddo

! Interpolate c,c`

        deallocate(mat_vfvf)
        allocate(mat_cfcc(xct%ncb_fi,xct%n2b_co,xct%nvb_fi,xct%nvb_fi))

        do iv=1,xct%nvb_fi
          do ivp=1,xct%nvb_fi
            mat_cfcc(1:xct%ncb_fi,1:xct%n2b_co,iv,ivp) = &
              MATMUL(conjg(dcck(1:xct%ncb_fi,1:xct%n2b_co,js_dcck)), &
              (mat_cccc(1:xct%n2b_co,1:xct%n2b_co,iv,ivp)))
          enddo
        enddo

        deallocate(mat_cccc)
        allocate(mat_cfcf(xct%ncb_fi,xct%ncb_fi,xct%nvb_fi,xct%nvb_fi))

        do iv=1,xct%nvb_fi
          do ivp=1,xct%nvb_fi
            mat_cfcf(1:xct%ncb_fi,1:xct%ncb_fi,iv,ivp) = &
              MATMUL((mat_cfcc(1:xct%ncb_fi,1:xct%n2b_co,iv,ivp)), &
              (dcckp(1:xct%n2b_co,1:xct%ncb_fi,js_dcckp)))
          enddo
        enddo

        deallocate(mat_cfcc)

! Reorder matrix

        do ivp=1,xct%nvb_fi
          do iv=1,xct%nvb_fi
            bse_fi(iv,1:xct%ncb_fi,ivp,1:xct%ncb_fi,js,jsp) = &
              mat_cfcf(1:xct%ncb_fi,1:xct%ncb_fi,iv,ivp)
          enddo
        enddo

        deallocate(mat_cfcf)

!      endif
    enddo
  enddo

  return
end subroutine interpolate

!> FHJ: Move a point (qq) to the 1st BZ (WS cell). Output vector is qq_bz, 
!! with length abs_qq2.
subroutine point_to_bz(crys, qq, qq_bz, abs_qq2)
  type(crystal), intent(in) :: crys
  real(DP), intent(in) :: qq(3)
  real(DP), intent(out) :: qq_bz(3)
  real(DP), intent(out) :: abs_qq2

  real(DP) :: abs_tmp, qq_tmp(3)
  integer :: i1, i2, i3

  !no push/pop, called too often
  abs_qq2 = INF
  do i1=-ncell,ncell+1
    qq_tmp(1) = qq(1) - i1
    do i2=-ncell,ncell+1
      qq_tmp(2) = qq(2) - i2
      do i3=-ncell,ncell+1
        qq_tmp(3) = qq(3) - i3
        abs_tmp = DOT_PRODUCT(qq_tmp, MATMUL(crys%bdot, qq_tmp))
        if (abs_tmp < abs_qq2) then
          abs_qq2 = abs_tmp
          qq_bz(:) = qq_tmp(:)
        endif
      enddo
    enddo
  enddo

end subroutine point_to_bz

end module intkernel_m
