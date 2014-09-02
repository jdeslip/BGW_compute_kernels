!===============================================================================
!
! Routines:
!
! (1) diag main         Originally by MLT       Last Edited: 5/12/2008 (JRD)
!
!     For more details see the README_absorption file.
!
!     Calculates the real and imaginary parts of the macroscopic dielectric
!     function starting from interaction matrix elements calculated by
!     the Kernel code. It uses interpolation in the matrix elements and
!     direct diagonalization of the Bethe-Salpeter equation.
!     Spin-polarized case implemented.
!
!     For more details, see:
!     Rohlfing & Louie, PRB 62:(8), 4927 (2000)
!     G. Strinati, Rivista del Nuovo Cimento 11:(12), 1 (1988)
!
!     Please report bugs to: jdeslip@civet.berkeley.edu
!     
!================================================================================

#include "f_defs.h"

subroutine diag(eqp,xct,flag,neig)

  use global_m
  use fullbz_m
  use genwf_m
  use intkernel_m
  use intwfn_m
  use misc_m
  implicit none

  type (eqpinfo), intent(inout) :: eqp
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(inout) :: flag
  integer, intent(inout) :: neig

  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (grid) :: kg_fi, kgq_fi
  type (kpoints) :: kp_fi, kpq_fi
  type (wavefunction) :: wfnc_fi
  type (wavefunction) :: wfnvq_fi
  type (work_genwf) :: work, workq
  type (int_wavefunction) :: intwfnc
  type (int_wavefunction) :: intwfnv

  character :: tmpstr*128,filename*20
  integer :: ii,jj,ncount,nmat,nblock,pblock
  integer :: ikb, icb, ivb
  integer :: ik,ikq,ikt,ikrq,ikcvs,ikcvsd,ic,iv,is
  integer :: iunit_c,iunit_v
  real(DP) :: vol,omega_plasma
  real(DP) :: tsec(2),tmin(2),tmax(2)
  
  character*16, allocatable :: routnam(:)
  integer, allocatable :: routsrt(:)
  integer, allocatable :: fi2co_wfn(:,:),indexq_fi(:)
  real(DP), allocatable :: egv(:), kco(:,:),cs(:),rdum(:)

  SCALAR, allocatable :: &
    dcc(:,:,:,:,:),dvv(:,:,:,:,:),s1(:),s1k(:,:,:),dummy(:), &
    egs(:,:),hqpcc(:,:,:),hqpvv(:,:,:),dipoles(:), &
    rdum2(:)
  !> (kcvs, k`c`v`s`), "A" block of BSE Hamiltonian
  SCALAR, allocatable :: hmtrx_A(:,:)
  !> (kcvs, k`c`v`s`), "B" block of BSE Hamiltonian, only if tda=.false.
  SCALAR, allocatable :: hmtrx_B(:,:)
  real(DP), allocatable :: intp_coefs(:,:)
  
  PUSH_SUB(diag)

! JRD: A Check for Finite Q

!      if (peinf%inode .eq. 0) then
!        write(6,*) 'nkpt_co = ', xct%nkpt_co
!      endif

  if(flag%vm.eq.2) then
    if (peinf%inode.eq.0) then
      write(0,*) 'WARNING: read_eps2_moments not supported in this diagonalization code. Ignoring keyword.'
    endif
    flag%vm=0
  endif
  
!--------------------------
! If eigenvalues.dat is available, read them and go straight to
! calculation of the absorption spectrum


  if (flag%spec.eq.1) then
    if (peinf%inode.eq.0) then
      call open_file(unit=14,file='eigenvalues.dat',form='formatted',status='old')
      read(14,'(10x,i8)') neig
      read(14,'(10x,e16.9)') vol
      read(14,'(10x,i8)') xct%nspin
      write(6,'(a,i8)') '# neig  = ', neig
      write(6,'(a,e16.9)') '# vol   = ', vol
      write(6,'(a,i8)') '# nspin = ', xct%nspin
      read(14,*)
      SAFE_ALLOCATE(cs, (neig))
      SAFE_ALLOCATE(egv, (neig))
      do ii=1,neig
        read(14,*) egv(ii),cs(ii)
      enddo
      call close_file(14)
    endif
    omega_plasma = 0.d0
    if (peinf%inode .eq. 0) then
      call absp(xct%eta,xct%sigma,xct%gamma,neig,cs,egv,vol,omega_plasma,flag)
    endif
    if (xct%iabsorp0 .eq. 0) then
      if (peinf%inode .eq. 0) write(6,*) 'Create absorption_eh.dat from eigenvalues.dat'
      if (peinf%inode .eq. 0) write(6,*) 'If you want absorption_noeh.dat add the noeh_only flag'
      call diag_end()
      POP_SUB(diag)
      return
    endif
  endif
!--------------------------
! Read wavefunctions on the fine grid

  call logit('Calling input')
  call timacc(2,1)
  call input(crys,gvec,kg_fi,kp_fi,syms,eqp,xct,flag, &
    omega_plasma,.true.,intwfnc)
  
! If there is no specified number of eigenvalues, calculate
! all eigenvalues

  nmat = xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin
  if (neig.eq.0) neig = nmat
  if ((neig.gt.nmat).or.(neig.le.0)) then
    write(tmpstr,'(a,2i6)') 'Incomprehensible request of eigenvalues : ', neig, nmat
    call die(tmpstr, only_root_writes = .true.)
  endif
  
  vol = xct%nktotal*crys%celvol
  if (peinf%inode.eq.0) then
    write(6,*) ' '
    write(6,*) 'More Job Parameters: '
    write(6,'(a,f32.14,a)') ' Crystal volume = ',vol,' a.u.'
    write(6,*) 'Number of valence bands = ',xct%nvb_fi
    write(6,*) 'Number of cond. bands   = ',xct%ncb_fi
    write(6,*) 'Number of spins   = ',xct%nspin
    write(6,'(a,f7.4,a)') ' Broadening: ',xct%eta,' eV'
    write(6,*) 'Number of eigenvalues to be computed = ',neig
    write(6,*) ' '
  endif
  call timacc(2,2)
      
  SAFE_ALLOCATE(indexq_fi, (xct%nkpt_fi))
  if (flag%vm.ne.1.or. .not. flag%read_dtmat) then
    call timacc(3,1)
    call logit('Calling input_q')
    call input_q(kp_fi,crys,gvec,kg_fi,kgq_fi,kpq_fi,syms,xct,indexq_fi,eqp,flag,intwfnv)
    call timacc(3,2)
  endif
      
! JRD: Don`t do this next section if only want absorption_noeh.dat

  if (xct%iabsorp0 .eq. 0) then

!------------------------------
! Calculate the transformation matrices from coarse grid wavefunctions
! FHJ: These are the final transformation coefs that will be used to interpolate
! the kernel. However, we might use an unrestricted version of dvv/dcc to
! interpolate eqp if xct%unrestricted_transf==.true..
    SAFE_ALLOCATE(dvv, (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
    SAFE_ALLOCATE(dcc, (xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
    SAFE_ALLOCATE(kco, (3,xct%nkpt_co))
    SAFE_ALLOCATE(fi2co_wfn, (xct%npts_intp_kernel,xct%nkpt_fi))
    SAFE_ALLOCATE(intp_coefs, (xct%npts_intp_kernel, xct%nkpt_fi))

    call logit('Calling intwfn')
    call timacc(4,1)
    call intwfn(kp_fi,crys,syms,xct,flag,gvec,kg_fi,kgq_fi,dcc,dvv,&
      kco,fi2co_wfn,indexq_fi,eqp,intwfnv,intwfnc,intp_coefs)
    call timacc(4,2)

  endif

  SAFE_DEALLOCATE_P(xct%ifmax)
  if (flag%vm.ne.1.or. .not. flag%read_dtmat) then
    ! otherwise, we did not call input_q to allocate it
    SAFE_DEALLOCATE_P(xct%ifmaxq)
  endif
  
!------------ Calculate the velocity (or momentum) matrix elements -------------

! Each PE calculates a small number of them. At the end, share data
!
! If flag%vm.eq.1, skip this part and just read the matrix elements
! from "vmtxel".
!
! nblock : size of cvs block in hmtrx


  call logit('Calculating v/p matrix elememts')
  nblock = xct%ncb_fi*xct%nvb_fi*xct%nspin
  if (xct%ipar .eq. 1) then
    peinf%nblockd=nblock
  else if (xct%ipar .eq. 2) then
    peinf%nblockd = xct%nvb_fi*xct%nspin
  else
    peinf%nblockd = xct%nspin
  endif
  nmat= xct%nkpt_fi*nblock
  SAFE_ALLOCATE(s1, (nmat))
  SAFE_ALLOCATE(s1k, (xct%ncb_fi,xct%nvb_fi,xct%nspin))
  s1= 0.d0
  s1k= 0.d0

  call timacc(10,1)

  if (flag%vm.eq.0) then
    do ikt=1, peinf%ikt(peinf%inode+1)
      ik = peinf%ik(peinf%inode+1,ikt)
      ikq = indexq_fi(ik)
      ikrq = kg_fi%indr(ik)
      
      call genwf(crys,gvec,kg_fi,syms,wfnc_fi,xct,ik,ik,work,intwfnc,is_cond = .true.)
      call genwf(crys,gvec,kgq_fi,syms,wfnvq_fi,xct,ik,ikq,workq,intwfnv,is_cond = .false.)

      if (flag%opr.eq.0) then
        call mtxel_v(wfnc_fi,wfnvq_fi,gvec,xct%qshift,wfnc_fi%nband,wfnvq_fi%nband,s1k)
      elseif (flag%opr.eq.1) then
        call mtxel_m(crys,wfnc_fi,wfnvq_fi,gvec,eqp,xct,wfnc_fi%nband,wfnvq_fi%nband,s1k,ik,.true.)
      endif
      do is=1,xct%nspin
        do ic=1,xct%ncb_fi
          do iv=1,xct%nvb_fi
            s1(bse_index(ik, ic, iv, is, xct)) = s1k(ic,iv,is)
          enddo
        enddo
      enddo
      SAFE_DEALLOCATE_P(wfnc_fi%cg)
      SAFE_DEALLOCATE_P(wfnc_fi%isort)
      SAFE_DEALLOCATE_P(wfnvq_fi%cg)
      SAFE_DEALLOCATE_P(wfnvq_fi%isort)
    enddo
    
    ! typedefs initializes all of these ikolds to 0
    if(work%ikold.ne.0) then
      SAFE_DEALLOCATE_P(work%cg)
      SAFE_DEALLOCATE_P(work%ph)
      SAFE_DEALLOCATE_P(work%ind)
      SAFE_DEALLOCATE_P(work%isort)
    endif
    if(workq%ikold.ne.0) then
      SAFE_DEALLOCATE_P(workq%cg)
      SAFE_DEALLOCATE_P(workq%ph)
      SAFE_DEALLOCATE_P(workq%ind)
      SAFE_DEALLOCATE_P(workq%isort)
    endif
! Share matrix elements

#ifdef MPI
    SAFE_ALLOCATE(dummy, (nmat))
    dummy = s1
    call MPI_ALLREDUCE(dummy,s1,nmat,MPI_SCALAR,MPI_SUM, &
      MPI_COMM_WORLD,mpierr)
    SAFE_DEALLOCATE(dummy)
#endif
    if (peinf%inode.eq.0) then
      write(6,*) ' writing matrix elements into vmtxel'
      call open_file(16,file='vmtxel',form='unformatted',status='replace')
      write(16) xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,flag%opr
      write(16) (s1(ikcvs),ikcvs=1,nmat)
      call close_file(16)

! If you want this file, you can get it by bsebinasc anyway.
!      call open_file(17,file='vmtxel.dat',form='formatted',status='replace')
!      write(17,*) xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,flag%opr
!      write(17,*) (s1(ikcvs),ikcvs=1,nmat)
!      call close_file(17)
    endif
    
  else ! ...if the matrix elements were calculated already
    
    if (peinf%inode.eq.0) then
      write(6,'(a)') ' reading matrix elements from vmtxel'
      call open_file(16,file='vmtxel',form='unformatted',status='old')
      read(16) ik,ic,iv,is,ii
      if (ik.ne.xct%nkpt_fi.or.ic.ne.xct%ncb_fi.or.iv.ne.xct%nvb_fi &
        .or.is.ne.xct%nspin.or.ii.ne.flag%opr) then
        write(0,'(a,5i6)') 'read  : ', ik,ic,iv,is,ii
        write(0,'(a,5i6)') 'needed: ', xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,flag%opr
        call die('parameter mismatch in vmtxel')
      endif
      
      read(16) (s1(ikcvs),ikcvs=1,nmat)
      call close_file(16)
    endif
    
#ifdef MPI
    call MPI_BCAST(s1,nmat,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
#endif

  endif

  call timacc(10,2)

  if (flag%vm.ne.1.or. .not. flag%read_dtmat) then
    call dealloc_grid(kgq_fi)
  endif

  if (xct%iwriteint.eq.1 .and. flag%vm == 0 .and. .not. flag%read_dtmat) then
    SAFE_DEALLOCATE_P(intwfnc%cgk)
    SAFE_DEALLOCATE_P(intwfnv%cgk)
    SAFE_DEALLOCATE_P(intwfnc%isort)
    SAFE_DEALLOCATE_P(intwfnv%isort)
  endif

! JRD: Now close the no longer needed wavefunction files

  if (xct%iwriteint.eq.0 .and. (flag%vm.ne.1 .or. .not. flag%read_dtmat)) then
    write(filename,'(a,i4.4)') 'INT_VWFNQ_', peinf%inode
    iunit_v = 128+(2*peinf%inode)+2
    call open_file(iunit_v, filename, status = 'old')
    call close_file(iunit_v, delete = .true.) ! files INT_VWFNQ_*
    
    write(filename,'(a,i4.4)') 'INT_CWFN_', peinf%inode
    iunit_c = 128+(2*peinf%inode)+1
    call open_file(iunit_c, filename, status = 'old')
    call close_file(iunit_c, delete = .true.) ! files INT_CWFN_*
  endif

! End Calculating Matrix Elements
!-------------------------------------------------------------------------------

!----------------------------
! Calculate the non-interacting spectrum. Only one PE works

  call timacc(9,1)
  call logit('Calling absp0')
  if (peinf%inode.eq.0) then
    call absp0(eqp,xct,s1,vol,omega_plasma,flag)
  end if
  call timacc(9,2)
  
  SAFE_DEALLOCATE_P(eqp%eclda)
  SAFE_DEALLOCATE_P(eqp%evlda)
  SAFE_DEALLOCATE(s1k)
  SAFE_DEALLOCATE(indexq_fi)

! JRD If we only want absorp0 - we finish here

  if (xct%iabsorp0 .eq. 1) then
    call diag_end()
    POP_SUB(diag)
    return
  endif

!------------ Build Hamiltonian Matrix --------------------------------------------

! Build Hamiltonian matrix. Diagonal part comes from the "non-interacting"
! quasiparticle Hamiltonians.  If the QP Greens function is diagonal,
! then these are just the quasiparticle energy differences on the
! diagonal.  The more general case is:
!
!       <cvk|H0|c'v'k'> = delta(k,k') *
!            [ <c|hqp|c'>*delta(v',v) - delta(c,c')*<v'|hqp|v> ]
!
! The rest of the Hamiltonian, which is the electron-hole interaction,
! comes from interpolation further below.

  call logit('Building non-interacting Hamiltonian')
  SAFE_ALLOCATE(hmtrx_A, (xct%nkpt_fi*nblock, peinf%nblocks*peinf%nblockd))
  hmtrx_A(:,:) = 0.0d0
  if (.not.xct%tda) then
    SAFE_ALLOCATE(hmtrx_B, (xct%nkpt_fi*nblock, peinf%nblocks*peinf%nblockd))
    hmtrx_B(:,:) = 0.0d0
  endif

! ikt loop. This loop is over proc owned k if ipar = 1, (k,c)  if 
! ipar = 2 and (k,c,v) if ipar = 3 

  do ikt=1,peinf%ibt(peinf%inode + 1)
    ik=peinf%ikb(peinf%inode+1,ikt)
    
    if (ik .eq. 0) then
      write(0,*) "Illegal value for ik",peinf%inode, ikt, ik
      call die("internal error in diag, ik = 0")
    endif
    
! Build <c|hqp|c'> and <v|hqp|v'> for this kpoint

    SAFE_ALLOCATE(hqpcc, (xct%ncb_fi,xct%ncb_fi,xct%nspin))
    SAFE_ALLOCATE(hqpvv, (xct%nvb_fi,xct%nvb_fi,xct%nspin))
    hqpcc = 0.0d0
    hqpvv = 0.0d0
    
! Put QP energies on diagonals of hqpcc and hqpvv to start

    do is=1,xct%nspin
      do ic=1,xct%ncb_fi
        hqpcc(ic,ic,is) = eqp%ecqp(ic,ik,is)
      enddo
      do iv=1,xct%nvb_fi
        hqpvv(iv,iv,is) = eqp%evqp(iv,ik,is)
      enddo
    enddo
    
! Read possible offdiagonal QP elements from "hqp.<ik>" file
! if it exists.  JRD: This is broken for now.  Someone should fix
! it in the future if they want to use it

    !if (ik.lt.10) then
    !  write(tmpstr,'(a,i1)') 'hqp.',ik
    !else if (ik.lt.100) then
    !  write(tmpstr,'(a,i2)') 'hqp.',ik
    !else if (ik.lt.1000) then
    !  write(tmpstr,'(a,i3)') 'hqp.',ik
    !else if (ik.lt.100000) then
    !  write(tmpstr,'(a,i5)') 'hqp.',ik
        !else
    !  write(0,*) 'too many kpoints for reading hqp'
    !endif
    !call open_file(9,file=tmpstr,form='formatted',status='old')
    !if (is.eq.0) then
    !  if (peinf%inode.eq.0) then
    !    write(6,*) 'Reading offdiagonal hqp from file ',tmpstr
    !    write(6,*) 'All values in eV'
    !  endif
    !  do
    !    read(9,*,end=999) nocc,ii,jj,x,y
    
    ! if ii and jj both refer to valence, states, put
    ! matrix element into hqpvv
    
    !    if ((ii<=nocc).and.(ii>nocc-xct%nvb_fi).and. &
    !    (jj<=nocc).and.(jj>nocc-xct%nvb_fi)) then
    !      if (peinf%inode.eq.0) write(6,'(a,2i5,2f20.10)') ' hqp(v,vp) = ',ii,jj,x,y
    !      ii=nocc-ii+1
    !      jj=nocc-jj+1
    !      is = 1
#ifdef CPLX
    !      hqpvv(ii,jj,is) = CMPLX(x,y)/ryd
#else
    !      hqpvv(ii,jj,is) = x/ryd
#endif
    !    else if ((ii>nocc).and.(ii<=nocc+xct%ncb_fi).and. &
    !    (jj>nocc).and.(jj<=nocc+xct%ncb_fi)) then
    !      if (peinf%inode.eq.0) write(6,'(a,2i5,2f20.10)') ' hqp(c,cp) = ',ii,jj,x,y
    !      ii=ii-nocc
    !      jj=jj-nocc
    !      is = 1
#ifdef CPLX
    !      hqpcc(ii,jj,is) = CMPLX(x,y)/ryd
#else
    !      hqpcc(ii,jj,is) = x/ryd
#endif
    !    endif
    !  enddo
    !999      call close_file(9)
    !  write(6,*)
    !endif ! if hqp.<ik> was found
    
    ! Now build hamiltonian from hqcc and hqvv
    
    do is=1,xct%nspin 
      do iv=1,xct%nvb_fi  
        do ic=1,xct%ncb_fi
          ikb=peinf%ikb(peinf%inode+1,ikt)
          if (xct%ipar .eq. 1) then
            ikcvsd = bse_index(ikt, ic, iv, is, xct)
            ivb=iv
            icb=ic
! JRD: if ipar is 2, the ikt loop above is over both k and c. So we don`t need
! the c loop here.
          else if (xct%ipar .eq. 2 .and. ic .eq. 1) then
            ikcvsd = bse_index(ikt, ic, iv, is, xct, ncband = 1)
            icb=peinf%icb(peinf%inode+1,ikt)
            ivb=iv
! JRD: if ipar is 3, the ikt loop above is over all over k,v and c. So we don`t need
! the c or v loops here.
          else if (xct%ipar .eq. 3 .and. iv .eq. 1 .and. ic .eq. 1) then
            ikcvsd = bse_index(ikt, ic, iv, is, xct, ncband = 1, nvband = 1)
            ivb=peinf%ivb(peinf%inode+1,ikt)
            icb=peinf%icb(peinf%inode+1,ikt)
          else
            cycle
          endif
          ikcvs = bse_index(ikb, icb, ivb, is, xct)
          
          hmtrx_A(ikcvs,ikcvsd) = hqpcc(icb,icb,is) - hqpvv(ivb,ivb,is)
        enddo
      enddo
    enddo
    SAFE_DEALLOCATE(hqpcc)
    SAFE_DEALLOCATE(hqpvv)
  enddo ! loop on k-points on this processor

!----------------------------
! Define the mapping of eigenvectors: the ii-th column of the matrix
! egs(:,:) stored in PE #ipe will contain the eigenvector of order
! peinf%peig(ipe,ii). The total number of eigenvectors stored in
! each processor is given by peinf%neig(1:peinf%npes).
! pblock >= maxval(peinf%neig(1:peinf%npes))


  pblock = neig/(peinf%npes*peinf%nblockd)
  if (pblock*peinf%npes*peinf%nblockd.lt.neig) pblock = pblock + 1
  pblock = pblock*peinf%nblockd
  SAFE_ALLOCATE(peinf%neig, (peinf%npes))
  SAFE_ALLOCATE(peinf%peig, (peinf%npes,pblock))
  peinf%neig=0
  peinf%peig=0
  ii=1
  do jj=1,neig
    if (ii.eq.peinf%npes+1) ii=1
    peinf%neig(ii)=peinf%neig(ii)+1
    peinf%peig(ii,peinf%neig(ii))=jj
    if (mod(jj,peinf%nblockd).eq.0) ii=ii+1
  enddo

!-----------------------------
! Interpolation scheme in the Kernel

  call logit('Calling intkernel')
  call timacc(5,1)
  if (xct%tda) then
    call intkernel(crys,kg_fi,kp_fi,syms,xct,hmtrx_A,dcc,dvv,kco,fi2co_wfn,flag,gvec,intp_coefs)
  else
    call intkernel(crys,kg_fi,kp_fi,syms,xct,hmtrx_A,dcc,dvv,kco,fi2co_wfn,flag,gvec,intp_coefs,hmtrx_B=hmtrx_B)
  endif
  call logit('Done intkernel')
  call timacc(5,2)
  SAFE_DEALLOCATE(fi2co_wfn)
  SAFE_DEALLOCATE(dcc)
  SAFE_DEALLOCATE(dvv)
  SAFE_DEALLOCATE(kco)

!#BEGIN_INTERNAL_ONLY
  !FHJ: This is for debugging purposes. Please keep this.
  if (.not.xct%tda.and.peinf%npes==1) then
    write(6,'(1x,a)') 'Dumping BSE Hamiltonian'
    call open_file(unit=666, file='hmtrx_A.dat', status='replace', form='formatted')
    call open_file(unit=667, file='hmtrx_B.dat', status='replace', form='formatted')
    write(666,'(2(i0,1x))') xct%nkpt_fi*nblock, peinf%nblocks*peinf%nblockd
    write(667,'(2(i0,1x))') xct%nkpt_fi*nblock, peinf%nblocks*peinf%nblockd
    do jj=1,peinf%nblocks*peinf%nblockd
#ifdef CMPLX
      write(666,'(2(es25.16,1x))') hmtrx_A(:, jj)
      write(667,'(2(es25.16,1x))') hmtrx_B(:, jj)
#else
      write(666,'(1(es25.16,1x))') hmtrx_A(:, jj)
      write(667,'(1(es25.16,1x))') hmtrx_B(:, jj)
#endif
    enddo
    call close_file(666)
    call close_file(667)
  endif
!#END_INTERNAL_ONLY

!--------------------------------
! Exact diagonalization

  call logit('Calling diagonalize')
  SAFE_ALLOCATE(egv, (neig))
  SAFE_ALLOCATE(egs, (nmat,pblock))
  call timacc(6,1)
  call diagonalize(pblock,neig,nmat,hmtrx_A,egv,egs)
  call timacc(6,2)
  
!--------------------------------
! Calculate transition matrix elements
! oscillator strength = 2 * cs / omega, as defined below

  call timacc(11,1)

  call logit('Computing transition matrix elements')
  SAFE_ALLOCATE(cs, (neig))
  SAFE_ALLOCATE(dipoles, (neig))
  cs = 0d0
  dipoles = ZERO

  do ii=1,peinf%neig(peinf%inode+1)
    jj= peinf%peig(peinf%inode+1,ii)
    dipoles(jj) = sum(MYCONJG( egs(1:nmat,ii) )*s1(1:nmat)) / sqrt(dble(xct%nspin))
    ! this factor is required to obtain the same transition matrix elements for the singlet excitons
    ! for nspin = 1 and nspin = 2, and is related to different normalization of the excitons
    ! see just after eq. (25) in Rohlfing and Louie, PRB 62, 4927 (2000)
    cs(jj) = ( abs(dipoles(jj)) )**2
  enddo

#ifdef MPI
  SAFE_ALLOCATE(rdum, (neig))
  SAFE_ALLOCATE(rdum2, (neig))
  rdum = cs
  rdum2 = dipoles
  call MPI_REDUCE(rdum,cs,neig,MPI_REAL_DP,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
  call MPI_REDUCE(rdum2,dipoles,neig,MPI_SCALAR,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
  SAFE_DEALLOCATE(rdum)
  SAFE_DEALLOCATE(rdum2)
#endif

  call timacc(11,2)

!-------------------------------
! Convert eigenvalues to eV and write them out

  egv(:) = egv(:)*ryd
  if (peinf%inode.eq.0) then
    if(any(egv(1:neig) < -TOL_Zero)) write(0,'(a)') "WARNING: There are negative excitation energies."
    call open_file(unit=14,file='eigenvalues.dat',form='formatted',status='replace')
    write(14,'(a,i8)') '# neig  = ', neig
    write(14,'(a,e16.9)') '# vol   = ', vol
    write(14,'(a,i8)') '# nspin = ', xct%nspin
    write(14,'(a)',advance='no') '#       eig (eV)   abs(dipole)^2'
    if (flag%krnl .ne. 0) then
      write(14,'(a)',advance='no') '      '
    else ! triplet transition matrix element = 0 if we consider spin overlap
      write(14,'(a)',advance='no') ' magn '
    endif
#ifdef CPLX
    write(14,'(a)') 'Re(dipole)      Im(dipole)'
#else
    write(14,'(a)') '    dipole'
#endif
    
    do ii=1,neig
      write(14,'(4e16.8)') egv(ii),cs(ii),dipoles(ii)
    enddo
    call close_file(14)
  endif
  
  SAFE_DEALLOCATE(s1)
  SAFE_DEALLOCATE(dipoles)
  
!------------------------------
! Calculate the absorption and density of excitonic states

  call timacc(12,1)

  call logit('Calling absp')
  if (peinf%inode.eq.0) then
    call absp(xct%eta,xct%sigma,xct%gamma,neig,cs,egv,vol,omega_plasma,flag)
  endif

  call timacc(12,2)

!------------------------------
! Write out eigenvectors to file if needed

  call timacc(13,1)

  if (flag%eig.ne.0) then
    call logit('Calling write_eigenvectors')
    call write_eigenvectors(xct,kg_fi,nblock,pblock, &
      neig,egv,egs,flag%eig)
  endif

  call timacc(13,2)
  
!-------------------------------
! Deallocate stuff

  SAFE_DEALLOCATE_P(peinf%neig)
  SAFE_DEALLOCATE_P(peinf%peig)
  call dealloc_grid(kg_fi)

  if (eqp%spl_tck%n>0) then
    SAFE_DEALLOCATE_P(eqp%spl_tck%t)
    SAFE_DEALLOCATE_P(eqp%spl_tck%c)
  endif

  SAFE_DEALLOCATE(cs)
  SAFE_DEALLOCATE(egv)
  SAFE_DEALLOCATE(hmtrx_A)
  if (.not.xct%tda) then
    SAFE_DEALLOCATE(hmtrx_B)
  endif
  SAFE_DEALLOCATE(egs)
  SAFE_DEALLOCATE_P(eqp%ecqp)
  SAFE_DEALLOCATE_P(eqp%evqp)
  SAFE_DEALLOCATE(intp_coefs)

  if (xct%iwritecoul .eq. 1 .and. peinf%inode .eq. 0) then
    call close_file(19) ! file vcoul
  endif

  call diag_end()

  POP_SUB(diag)
  return

contains

  subroutine diag_end()

    PUSH_SUB(diag.diag_end)

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

!--------------------------------
! Time accounting

    SAFE_ALLOCATE(routnam, (60))
    routnam(1)='TOTAL:'
    routnam(2)='INPUT:'
    routnam(3)='INPUT_Q:'
    routnam(4)='INTWFN:'
    routnam(5)='INTKERNEL:'
    routnam(6)='DIAGONALIZE:'
    routnam(7)='EPSDIAG:'
    routnam(8)='EPS COMM:'
    routnam(9)='ABSORP0:'
    routnam(10)='VMTXEL:'
    routnam(11)='TRANS MTXEL:'
    routnam(12)='ABSORP:'
    routnam(13)='WRITE EIG:'
    routnam(41)='IW Input_co:'
    routnam(42)='IW Interp:'
    routnam(43)='IW Genwf:'
    routnam(44)='IW Gwnwf_Co:'
    routnam(45)='IW Mtxel_t:'
    routnam(46)='IW Write:'
    routnam(47)='IW Reduce:'
    routnam(51)='IK Setup:'
    routnam(52)='IK C-Check:'
    routnam(53)='IK Input:'
    routnam(54)='IK Inteps:'
    routnam(55)='IK Vcoul:'
    routnam(56)='IK Cache:'
    routnam(57)='IK Interp:'
    routnam(58)='IK Sum:'
    SAFE_ALLOCATE(routsrt, (28))
    routsrt=(/ (ii,ii=2,13), (ii,ii=41,47), (ii,ii=51,58), 1 /)
  
    call timacc(1,2)
    if(peinf%inode.eq.0) then
      write(6,*)
      write(6,9000) 'CPU [s]','WALL [s]','#'
      write(6,*)
    endif
    
    do ii=1,ubound(routsrt, 1)
      call timacc(routsrt(ii),3,tsec,ncount)
#ifdef MPI
      call MPI_ALLREDUCE(tsec,tmin,2,MPI_REAL_DP,MPI_MIN,MPI_COMM_WORLD,mpierr)
      call MPI_ALLREDUCE(tsec,tmax,2,MPI_REAL_DP,MPI_MAX,MPI_COMM_WORLD,mpierr)
#else
      tmin = tsec
      tmax = tsec
#endif
      if(peinf%inode==0) then
        if (ii>1) then
          if (routsrt(ii)-routsrt(ii-1)/=1) write(6,*)
        endif
        write(6,9001) routnam(routsrt(ii)),tmin(1),tmin(2),ncount
        write(6,9002) tsec(1),tsec(2)
        write(6,9003) tmax(1),tmax(2)
      endif
    enddo
    
9000 format(22x,a13,  3x,a13,  3x,a8)
9001 format(1x,a16,'(min.)',f13.3,3x,f13.3,3x,i8)
9002 format(  17x,'(PE 0)',f13.3,3x,f13.3)
9003 format(  17x,'(max.)',f13.3,3x,f13.3)

    POP_SUB(diag.diag_end)
    return
  end subroutine diag_end

end subroutine diag
