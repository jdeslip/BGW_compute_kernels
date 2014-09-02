program kernel_BSE

  use common_m
  use intkernel_m
  implicit none

  type(crystal) :: crys
  type(grid) :: kg
  type (xctinfo) :: xct
  complex(DPC), allocatable :: dcc(:,:,:,:,:), dvv(:,:,:,:,:)
  complex(DPC), allocatable :: hmtrx(:,:)
  complex(DPC) :: tmp(2,2)
  real(DP), allocatable :: intp_coefs(:,:)
  integer, allocatable :: fi2co_wfn(:,:)

  integer :: nmat, kgrid_co(3), kgrid_fi(3), ik, is, iv, ic, ikcvs
  integer :: i1, i2, i3, ii

  ! FHJ - Kernel initialization
  xct%npts_intp_kernel = 1
  kgrid_co(1) = 2
  kgrid_co(2) = 2
  kgrid_co(3) = 2
  xct%nkpt_co = product(kgrid_co)
  kgrid_fi(1) = 4
  kgrid_fi(2) = 4
  kgrid_fi(3) = 4
  xct%nkpt_fi = product(kgrid_fi)
  xct%nktotal = xct%nkpt_fi
  xct%nspin = 1
  xct%nvb_fi = 2
  xct%ncb_fi = 2
  xct%nvb_co = 4
  xct%ncb_co = 4
  xct%n1b_co = xct%nvb_co
  xct%n2b_co = xct%ncb_co
  xct%idimensions = 3
  xct%skipinterp = .false.

  write(6,*)
  write(6,*) 'Initializing default params'

  crys%celvol = 1d0
  crys%bdot = 0d0
  crys%bdot(1,1) = 1d0
  crys%bdot(2,2) = 1d0
  crys%bdot(3,3) = 1d0
  xct%wcoul0 = 10d0

  nmat = xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin
  allocate(hmtrx(nmat,nmat))
  hmtrx(:,:) = (0d0, 0d0)
  do ik=1, xct%nkpt_fi
    do is=1,xct%nspin 
      do iv=1,xct%nvb_fi  
        do ic=1,xct%ncb_fi
          ikcvs = bse_index(ik, ic, iv, is, xct)
          hmtrx(ikcvs,ikcvs) = 1d0 + ic*1d0 + iv*1d0 + (ik-xct%nkpt_fi*0.5d0)*1d-3
        enddo
      enddo
    enddo
  enddo ! loop on k-points on this processor
  hmtrx(:,:) = (0d0, 0d0)

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

  write(6,*)
  write(6,*) 'Calling intkernel'
  call intkernel(crys,kg,xct,hmtrx,dcc,dvv,fi2co_wfn,intp_coefs)

  write(6,*)
  write(6,*) 'Summary of resulting matrix:'
  write(6,*) 'Sum(H) = ', sum(hmtrx)
  write(6,*) 'Sum(|real(H)|) = ', sum(abs(dble(hmtrx)))
  write(6,*) 'Sum(|imag(H)|) = ', sum(abs(aimag(hmtrx)))
  write(6,*) 'Norm2(real(H)+imag(H)) = ', norm2(dble(hmtrx)+aimag(hmtrx))
  write(6,*) 'H(1:2,1:2) = '
  tmp = hmtrx(1:2, 1:2)
  write(6,*) tmp
  write(6,*)

end program kernel_BSE
