!==============================================================================
! kernel_chi_sum.x : Kernel for chi summation in BGW.
! By Derek Vigil-Fowler, with credit to Felipe da Jornada at many steps.
!
!==============================================================================

program chi_summation

  use blas_m
  use timing_m

  implicit none

    integer :: icurr,ntot,ntot2,itot,ifreq_para
    integer :: ipe, ilimit, jj, iv, irk, it, ispin,im,mytot
    integer :: i_myband,grp_mtxel_start,tmp_iv,im_proc
    complex(kind((1.0d0,1.0d0))) :: negfact,tmp(2,2)
    real(kind(1.0d0)) :: zvalue
    integer, allocatable :: tmprowindex(:),tmpcolindex(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: tmprowph(:),tmpcolph(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: gmeRDyntempr(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: edenDRtemp(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: edenDAtemp(:)
    real(kind(1.0d0)) :: fact
    real(kind(1.0d0)), parameter :: ryd = 1.360569253d1
    real(kind(1.0d0)), parameter :: Tol_Zero = 1.0d-12
    ! DVF : variables replacing those in structures in standard BGW
    ! Also some new variables introduced for initializing arrays
    ! that are normally provided as input to chi_summation.
    integer :: nfreq,nrk,nval,ncond,os_para_freqs,os_nfreq_para,&
               nspin,ifreq,nthreads,ic,npes,inode,tid,ig,ik,&
               OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,narg,nprc_nmtx,&
               ncount,ii
    real(kind(1.0d0)) :: cond_en, val_en,tsec(2)
    integer, allocatable :: nst(:),indt(:,:,:)
    complex(kind((1.0d0,1.0d0))), allocatable :: pht(:,:,:),gme(:,:,:,:,:,:),dFreqBrd(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: chiRDyn(:,:,:,:),chiADyn(:,:,:,:)
    real(kind((1.0d0,1.0d0))), allocatable :: edenDyn(:,:,:,:,:),dFreqGrid(:)
    character*16, allocatable :: routnam(:)
    character*16 :: arg
    integer, allocatable :: routsrt(:)
    !> DWV: below are some temporary buffers needed for the chi summation. They are
    !! described in detail in this comment.
    !! gme = g-vector matrix element
    !! gmetempX where X = n,r,c are temporary arrays for static calculations
    !! n = normalized by the proper factor used in BGW
    !! r = row, meaning the matrix elements for all bands (nv*nc*nk) that the proc owns
    !! for the G-vectors in the rows of chi currently being computed
    !! c = column, the matrix elements for all bands (nv*nc*nk) that the proc owns
    !! for the G`-vectors in the rows of chi currently being computed
    !! the RDyn arrays are needed for full-frequency (FF) calculations, real and complex
    !! while the Adyn arrays are needed only for complex FF calculations
    !! r2 is used in FF with matrix communication because you have a different denominator for 
    !! each frequency. Only the r2 array (not the r) array is used for element communication
    !! the denominators are built into the gme`s for static calculations
    !! eden arrays hold the energy denominators for FF
    !! chilocal holds the contribution of a given processor to the GG` chunk of epsilon
    !! being computed
    !! chilocal2 holds the MPI reduced GG` chunk of epsilon being computed 
    complex(kind((1.0d0,1.0d0))), allocatable :: chilocalRDyn(:,:,:)
    complex(kind((1.0d0,1.0d0))), allocatable :: chilocal2RDyn(:,:,:,:)
    complex(kind((1.0d0,1.0d0))), allocatable :: chiRDyntmp(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: gmeRDyntempn(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: gmeRDyntempr2(:,:,:)
    complex(kind((1.0d0,1.0d0))), allocatable :: gmeRDyntempc(:,:)

    complex(kind((1.0d0,1.0d0))), allocatable :: chiADyntmp(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: chilocalADyn(:,:,:)
    complex(kind((1.0d0,1.0d0))), allocatable :: chilocal2ADyn(:,:,:,:)
    complex(kind((1.0d0,1.0d0))), allocatable :: gmeADyntempn(:)
    complex(kind((1.0d0,1.0d0))), allocatable :: gmeADyntempr2(:,:,:)
    complex(kind((1.0d0,1.0d0))), allocatable :: gmeADyntempc(:,:) 

    ! Here's where we read in and assign all the variables that will be 
    ! needed for the computation. Borrowed heavily from BSE kernel of 
    ! FHJ. 
    write(6,*)
    write(6,*) 'BerkeleyGW - Kernel for the chi summation'
    write(6,*) '==========================================================='
    write(6,*)
    narg = iargc()
    if (narg/=5) then
      write(0,*) 'Usage:'
      write(0,*) ' kernel_chi_sum.x nv nc ng nk nfreq'
      write(0,*)
      write(0,*) 'Short example (<1s):'
      write(0,*) './kernel_chi_sum.x  4 100 100 5 15'
      write(0,*)
      write(0,*) 'Long example (<30s):'
      write(0,*) './kernel_chi_sum.x  4 400 100 30 20'
      write(0,*)
      stop
    endif

    ! DVF - Kernel initialization
    ! ---------------------------------------------------------------------------

    call getarg(1, arg)
    read(arg,*) nval
    call getarg(2, arg)
    read(arg,*) ncond
    call getarg(3, arg)
    read(arg,*) nprc_nmtx
    call getarg(4, arg)
    read(arg,*) nrk
    call getarg(5, arg)
    read(arg,*) nfreq

    write(6,*) 'Input parameters'
    write(6,*) '----------------'
    write(6,666) 'number of valence bands', nval
    write(6,666) 'number of conduction bands', ncond
    write(6,666) 'number of gvectors in polarizability', nprc_nmtx
    write(6,666) 'number of kpoints', nrk
    write(6,666) 'number of frequencies', nfreq

!$omp parallel private(tid) shared(nthreads)
    tid = OMP_GET_THREAD_NUM()
    if (tid==0) then
      nthreads = OMP_GET_NUM_THREADS()
      write(6,666) 'Number of threads', nthreads
    endif
!$omp end parallel
    666 format(1x,a,' = ',i0)
    write(6,*)

    call timacc(0,0)
    call timacc(1,1)

    write(6,*) 'Initializing default params'
    call timacc(2,1)

    os_para_freqs=1
    os_nfreq_para=nfreq
    nspin=1
    npes=1
    inode=1
    allocate(nst(nrk))
    nst(:)=1

    allocate(dFreqGrid(nfreq))
    allocate(dFreqBrd(nfreq))

! Frequency grid and broadening. This is reasonable for the
! semiconductor model I`m using for the bandstructure

    do ifreq=1,nfreq
      dFreqGrid(ifreq)=real((ifreq-1))/real(nfreq)*100.0
      dFreqBrd(ifreq)= cmplx(0d0,5d-1)
    enddo

! We`ll assume we have some 1D parabolic bands here, with
! the first kpoint being k=0 and then we go out farther 
! in the BZ as the kpoint index increases. This BS is 
! not physical, but it`s easy to contrive and shouldn't
! cause anything to blow up

    allocate(edenDyn(nval,ncond,nspin,nrk,os_para_freqs))    

    do iv=1,nval
      do ic=1,ncond
        do ik=1,nrk
          val_en  = (iv-1)*1d0 + real(ik-1)**2/real(nrk-1)**2*1d1
          cond_en = nval*1.1d1 + 2d0 + ic*5d-1 + real(ik-1)**2/real(nrk-1)**2*5d0
          edenDyn(iv,ic,1,ik,1)= (cond_en-val_en)/ryd
        enddo
      enddo
    enddo 

! Matrix elements are not easy to model, so I won`t really try. 
! We`ll assume this is the small q-point, so the states are close
! to orthogonal. We'll also assume weak k-point dependence and that
! the higher conduction states are more orthogonal to the valence 
! states than the lower ones.

    allocate(gme(nprc_nmtx,ncond,nval,nspin,nrk,os_para_freqs))

    do ik=1,nrk
      do iv=1,nval
        do ic=1,ncond
          do ig=1,nprc_nmtx
            if( (ig-ic-iv) .eq. 0 .or. (ic-ig-iv) .eq. 0) then
              gme(ig,ic,iv,1,ik,1)=1d-1
            else
              gme(ig,ic,iv,1,ik,1)= 1d-2/( abs( abs( real(ig)-real(ic) )-real(iv) )*2d-1 )
            endif
          enddo
        enddo
      enddo
    enddo

! We're pretending there`s no symmetry, as seen by the fact that
! we put 1 k-point in the star of all k-points. So, all phases
! are 1 and the index array is just the identity map. 
    
    allocate(indt(nprc_nmtx,nst(1),nrk))
    allocate(pht(nprc_nmtx,nst(1),nrk))

    do ig=1,nprc_nmtx
      do ik=1,nrk
        indt(ig,1,ik)=ig
        pht(ig,1,ik) = cmplx(1d0,0d0)
      enddo
    enddo

    call timacc(2,2)

    call timacc(3,1)
    call timacc(6,1)
    allocate(chilocal2RDyn(nprc_nmtx,nprc_nmtx,os_nfreq_para,1))
    chilocal2RDyn=0
    allocate(chilocal2ADyn(nprc_nmtx,nprc_nmtx,os_nfreq_para,1))
    chilocal2ADyn=0

    allocate(gmeRDyntempn(nprc_nmtx))
    allocate(chiRDyntmp(os_nfreq_para))
    allocate(gmeADyntempn(nprc_nmtx))
    allocate(chiADyntmp(os_nfreq_para))
    call timacc(6,2)

    ntot=0
    ntot2=0

    ntot = nval*ncond
    do irk = 1, nrk
      ntot2=ntot2 + nst(irk)
    enddo
    ntot=ntot*ntot2

    fact= 1d-3 ! This is a reasonable value for this number 

    negfact = -1D0*fact

    grp_mtxel_start=inode !for non-parallel freq case

    do ipe = 1, npes

      call timacc(6,1)
      allocate(gmeRDyntempr2(nprc_nmtx,ntot,os_nfreq_para))
      allocate(gmeRDyntempc(ntot,nprc_nmtx))
      allocate(chilocalRDyn(nprc_nmtx,nprc_nmtx,os_nfreq_para))
      chilocalRDyn=0
      allocate(gmeADyntempr2(nprc_nmtx,ntot,os_nfreq_para))
      allocate(gmeADyntempc(ntot,nprc_nmtx))
      allocate(chilocalADyn(nprc_nmtx,nprc_nmtx,os_nfreq_para))
      chilocalADyn=0
      call timacc(6,2)

      do ispin = 1 , nspin

        itot = 0

        call timacc(4,1)

        allocate(tmprowindex(nprc_nmtx))
        allocate(tmpcolindex(nprc_nmtx))
        allocate(tmprowph(nprc_nmtx))
        allocate(tmpcolph(nprc_nmtx))

        do im=1,os_para_freqs                              ! im labels which member of the mtxel comm group are you
          im_proc= 1                                       ! im_proc gives this mtxel comm group member`s global
          do irk = 1, nrk                                  ! proc number
            do it = 1, nst(irk)

              do icurr=1,nprc_nmtx
                tmprowindex(icurr) = indt(icurr,it,irk)
                tmprowph(icurr) = pht(icurr,it,irk)
              enddo
              do icurr=1,nprc_nmtx
                tmpcolindex(icurr) = indt(icurr,it,irk)
                tmpcolph(icurr) = pht(icurr,it,irk)
              enddo
              do iv = 1,nval
                tmp_iv = iv

                if (.true.) then
                  ilimit = ncond
                else
                  ilimit = 0
                endif

                !$OMP PARALLEL private (mytot,zvalue,edenDRtemp,edenDAtemp,gmeRDyntempr,jj,icurr)
                allocate(edenDAtemp(os_nfreq_para))

                allocate(gmeRDyntempr(nprc_nmtx))
                allocate(edenDRtemp(os_nfreq_para))
                !$OMP DO
                do i_myband = 1, ilimit
                  mytot = itot + i_myband
                  zvalue = edenDyn(iv,i_myband,ispin,irk,im)
                  if(.true.) then
                    !this is when the lin_denominator mode is not active.
                    do jj=1,nfreq,os_para_freqs
                      ifreq_para=(jj+os_para_freqs-1)/os_para_freqs
                      if (abs(zvalue) .gt. Tol_Zero) then
                        edenDRtemp(ifreq_para)= -0.5d0*( &
                          1d0/(zvalue-(dFreqBrd(jj)+dFreqGrid(jj))/ryd)+ &
                          1d0/(zvalue+(dFreqBrd(jj)+dFreqGrid(jj))/ryd))
                      else
                        edenDRtemp(ifreq_para)= 0D0
                      endif
                    enddo

                  endif 

                  if(.true.) then
                    !this is when the lin_denominator mode is not active.
                    do jj=1,nfreq,os_para_freqs
                      ifreq_para=(jj+os_para_freqs-1)/os_para_freqs
                      if (abs(zvalue) .gt. Tol_Zero) then
                        edenDAtemp(ifreq_para)= -0.5d0*( &
                          1d0/(zvalue-(-dFreqBrd(jj)+dFreqGrid(jj))/ryd)+ &
                          1d0/(zvalue+(-dFreqBrd(jj)+dFreqGrid(jj))/ryd))
                      else
                        edenDAtemp(ifreq_para)= 0D0
                      endif
                    enddo

                  endif

                  do icurr=1,nprc_nmtx
                    gmeRDyntempr(icurr)=gme(tmprowindex(icurr), &
                      i_myband,tmp_iv,ispin,irk,im) * tmprowph(icurr)
                  enddo
                  do jj = 1,nfreq,os_para_freqs
                    ifreq_para=(jj+os_para_freqs-1)/os_para_freqs
                    gmeRDyntempr2(:,mytot,ifreq_para)=gmeRDyntempr(:)*edenDRtemp(ifreq_para)
                  enddo
                  do icurr=1,nprc_nmtx
                    gmeRDyntempc(mytot,icurr) = &
                      conjg( gme(tmpcolindex(icurr),i_myband,tmp_iv,ispin,irk,im) * tmpcolph(icurr) )
                  enddo
                  do jj = 1,nfreq,os_para_freqs
                    ifreq_para=(jj+os_para_freqs-1)/os_para_freqs
! POSSIBLY BAD FOR OPENMP
                    gmeADyntempr2(:,mytot,ifreq_para)=gmeRDyntempr(:)*edenDAtemp(ifreq_para)
                  enddo
                  do icurr=1,nprc_nmtx
! PROBABLY BAD FOR OPENMP
                    gmeADyntempc(mytot,icurr)= &
                      conjg( gme(tmpcolindex(icurr),i_myband,tmp_iv,ispin,irk,im) * tmpcolph(icurr) )
                  enddo
                enddo ! i_myband
                !$OMP END DO
                deallocate(edenDRtemp)
                deallocate(edenDAtemp)
                deallocate(gmeRDyntempr)
                !$OMP END PARALLEL
                itot = itot+ilimit
              enddo ! iv
            enddo ! it
          enddo ! irk
        enddo ! im

        deallocate(tmprowindex)
        deallocate(tmpcolindex)
        deallocate(tmprowph)
        deallocate(tmpcolph)

        call timacc(4,2)          

        !Do the zgemm`s
        if(ntot > 0) then
          do jj =1,nfreq,os_para_freqs
            ifreq_para=(jj+os_para_freqs-1)/os_para_freqs
            call timacc(5,1)
            call zgemm('n','n',nprc_nmtx,nprc_nmtx,ntot, &
              negfact,gmeRDyntempr2(:,:,ifreq_para),nprc_nmtx,gmeRDyntempc(:,:),ntot,&
              (0D0,0D0),chilocalRDyn(:,:,ifreq_para),nprc_nmtx)
            call timacc(5,2)
            call timacc(5,1)
            call zgemm('n','n',nprc_nmtx,nprc_nmtx,ntot, &
              negfact,gmeADyntempr2(:,:,ifreq_para),nprc_nmtx,gmeADyntempc(:,:),ntot,&
              (0D0,0D0),chilocalADyn(:,:,ifreq_para),nprc_nmtx)
            call timacc(5,2)
          enddo
        endif

        call timacc(8,1)

        chilocal2RDyn(:,:,:,ispin)=chilocalRDyn(:,:,:)
        chilocal2ADyn(:,:,:,ispin)=chilocalADyn(:,:,:)

        call timacc(8,2)

      enddo ! ispin
      call timacc(6,1)
      deallocate(chilocalRDyn)
      deallocate(gmeRDyntempr2)
      deallocate(gmeRDyntempc)
      deallocate(chilocalADyn)
      deallocate(gmeADyntempr2)
      deallocate(gmeADyntempc)
      call timacc(6,2)
    enddo ! ipe

    call timacc(8,1)

    allocate(chiRDyn(nfreq,nprc_nmtx,nprc_nmtx,1))
    allocate(chiADyn(nfreq,nprc_nmtx,nprc_nmtx,1))

    do ispin =1, nspin
      do jj=1,nfreq,os_para_freqs
        ifreq_para=(jj+os_para_freqs-1)/os_para_freqs
        chiRDyn(ifreq_para,:,:,ispin) = chilocal2RDyn(:,:,ifreq_para,ispin)
        chiADyn(ifreq_para,:,:,ispin) = chilocal2ADyn(:,:,ifreq_para,ispin)
      enddo ! jj
    enddo ! ispin

    call timacc(8,2)

    call timacc(6,1)
    deallocate(chilocal2RDyn)
    deallocate(chilocal2ADyn)

    deallocate(gmeRDyntempn)
    deallocate(chiRDyntmp)
    deallocate(gmeADyntempn)
    deallocate(chiADyntmp)
    call timacc(6,2)
    call timacc(3,2)

    deallocate(nst)
    deallocate(dFreqGrid)
    deallocate(dFreqBrd)
    deallocate(edenDyn)    
    deallocate(gme)
    deallocate(indt)
    deallocate(pht)


    ! DVF - Summarize results (borrowed from FHJ again)
    ! ---------------------------------------------------------------------------

    call timacc(7,1)
    write(6,*)
    write(6,*) 'Summary of resulting matrix:'
    write(6,*) '----------------------------'
    write(6,*) 'Sum(chi0R(omega=0)) = ', sum(chiRDyn(1,:,:,1))
    write(6,*) 'Sum(|real(chi0R(omega=0))|) = ', sum(abs(dble(chiRDyn(1,:,:,1))))
    write(6,*) 'Sum(|imag(chi0R(omega=0))|) = ', sum(abs(aimag(chiRDyn(1,:,:,1))))
    write(6,*) 'Norm2(real(chi0R(omega=0))+imag(chi0R(omega=0))) = ',norm2(dble(chiRDyn(1,:,:,1))+aimag(chiRDyn(1,:,:,1)))
    write(6,*) 'chi0R(1:2,1:2)(omega=0) = '
    tmp = chiRDyn(1,1:2,1:2,1)
    write(6,*) tmp
    write(6,*)
    call timacc(7,2)

    deallocate(chiRDyn)
    deallocate(chiADyn)

    ! DVF - Timing (ditto)
    ! ---------------------------------------------------------------------------

    allocate(routnam(8))
    routnam(1)='TOTAL:'
    routnam(2)='SETUP:'
    routnam(3)='Chi Sum'
    routnam(4)='CS Prep'
    routnam(5)='CS ZGEMM'
    routnam(6)='CS Alloc'
    routnam(7)='SUMMARIZE:'
    routnam(8)='CS Array Copy'
    allocate(routsrt(8))
    routsrt=(/ 2, 3, 4,5, 6, 8, 7, 1 /)

    call timacc(1,2)
    write(6,*)
    write(6,9000) 'CPU (s)','WALL (s)','#'
    write(6,*)
    do ii=1,ubound(routsrt, 1)
      call timacc(routsrt(ii),3,tsec,ncount)
      if (ii>1) then
        if (abs(routsrt(ii)-routsrt(ii-1))>10) write(6,*)
      endif
      write(6,9001) routnam(routsrt(ii)),tsec(1)/nthreads,tsec(2),ncount
    enddo
    write(6,*)

    deallocate(routnam)
    deallocate(routsrt)

9000 format(22x,a13,  3x,a13,  3x,a8)
9001 format(1x,a16,5x,f13.3,3x,f13.3,3x,i8)

end program chi_summation 
