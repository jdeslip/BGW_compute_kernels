!=============================================================================================
!
! This Kernel Represents the Full-Frequency Self-Energy Summations in BerkeleyGW
!
! Run like:
!
! ffkernel.x <number_bands> <number_valence_bands> <number_plane_waves> <number_plane_waves_per_mpi_task> <number_frequencies> <number_evalutation_energies>
!
! Example of bandwidth bound run (large system fewer frequencies):
!
!  ffkernel.x 8 2 10000 600 400 200 
!
! Example of compute bound run (small system asking for many frequencies):
!  ffkernel.x 8 2 1000 60 20000 30000
!
!==============================================================================================

program ffkernel

implicit none

integer :: NTHREADS,TID,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
integer :: n1
integer :: ispin, iw, ifreq, ijk, iwlda
integer :: number_bands,nvband,ncouls,ngpown,nfreqeval,nFreq
integer :: my_igp, indigp, ig, igp, igmax
integer :: ggpsum
integer, allocatable :: indinv(:), inv_igp_index(:)
complex(kind((1.0d0,1.0d0))) :: achstemp,achxtemp,matngmatmgp,matngpmatmg,mygpvar1,mygpvar2,schstemp,schs,sch,ssx,ssxt,scht
complex(kind((1.0d0,1.0d0))), allocatable :: aqsmtemp(:,:), aqsntemp(:,:), I_epsR_array(:,:,:), I_epsA_array(:,:,:),matngmatmgpD(:,:),matngpmatmgD(:,:),dFreqBrd(:)
complex(kind((1.0d0,1.0d0))) :: schD,achsDtemp,schsDtemp
complex(kind((1.0d0,1.0d0))), allocatable :: asxDtemp(:),achDtemp(:),ach2Dtemp(:),achDtemp_cor(:),achDtemp_corb(:)
complex(kind((1.0d0,1.0d0))), allocatable :: schDi(:),schDi_cor(:),schDi_corb(:),sch2Di(:),ssxDi(:)
complex(kind((1.0d0,1.0d0))) :: ssxDit,ssxDitt,ssxDittt,schDt,schDtt,sch2dt,sch2dtt,I_epsRggp_int, I_epsAggp_int
complex(kind((1.0d0,1.0d0))), allocatable :: schDt_array(:),schDt_matrix(:,:)
complex(kind((1.0d0,1.0d0))) :: schDttt,schDttt_cor
complex(kind((1.0d0,1.0d0))) :: schDt_avg, schDt_right, schDt_left, schDt_lin, schDt_lin2, schDt_lin3
complex(kind((1.0d0,1.0d0))) :: cedifft_coh,cedifft_cor
real(kind(1.0d0)) :: cedifft_zb,intfact,cedifft_zb_left,cedifft_zb_right
real(kind(1.0d0)) :: e_n1kq, e_lk, dw, pref_zb
real(kind(1.0d0)) :: tol,fact1,fact2
real(kind(1.0d0)) :: wx,occ,occfact
real(kind(1.0d0)), allocatable :: ekq(:,:), vcoul(:), wxi(:), dFreqGrid(:), pref(:)
real(kind(1.0d0)) :: limitone,limittwo
real(kind(1.0d0)) :: starttime, endtime
real(kind(1.0d0)) :: starttime_stat, endtime_stat
real(kind(1.0d0)) :: time_stat
real(kind(1.0d0)) :: starttime_sx, endtime_sx
real(kind(1.0d0)) :: time_sx 
real(kind(1.0d0)) :: starttime_ch, endtime_ch
real(kind(1.0d0)) :: time_cha,time_chb,time_chc 
real(kind(1.0d0)) :: freqevalmin, freqevalstep
logical :: flag_occ

CHARACTER(len=32) :: arg

time_stat = 0D0
time_sx = 0D0
time_cha = 0D0
time_chb = 0D0
time_chc = 0D0

!$OMP PARALLEL PRIVATE(NTHREADS, TID)
      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        NTHREADS = OMP_GET_NUM_THREADS()
        PRINT *, 'Number of threads = ', NTHREADS
      END IF
!$OMP END PARALLEL

! We start off in the body of loop over the various tasks. Each MPI task has communicated data it owns to everyone

! These should be read in

      CALL getarg(1, arg)
      READ(arg,*) number_bands
      CALL getarg(2, arg)
      READ(arg,*) nvband
      CALL getarg(3, arg)
      READ(arg,*) ncouls
      CALL getarg(4, arg)
      READ(arg,*) ngpown
      CALL getarg(5, arg)
      READ(arg,*) nFreq
      CALL getarg(6, arg)
      READ(arg,*) nfreqeval
      
      write(6,*)"number_bands = ",number_bands
      write(6,*)"nvband       = ",nvband
      write(6,*)"ncouls       = ",ncouls
      write(6,*)"ngpown       = ",ngpown  ! this is also the mpi level parallel
      write(6,*)"nFreq        = ",nFreq
      write(6,*)"nfreqeval    = ",nfreqeval

! ngpown = ncouls / (number of mpi tasks)

      freqevalmin = 0D0
      freqevalstep = 0.5D0

      ggpsum = 2
      e_lk = 10D0
      dw = 1D0
      pref_zb = 0.5D0 / 3.14D0

      ALLOCATE(vcoul(ncouls))
      vcoul = 1D0

      ALLOCATE(ekq(number_bands,1))
      dw = -10D0
      do ijk = 1, number_bands
        ekq(ijk,1) = dw 
        dw = dw + 1D0
      enddo

      write(6,*)"allocating aqs*temp, size = ",2*ncouls*number_bands*16/1024," kbytes"
      ALLOCATE(aqsntemp(ncouls,number_bands))
      ALLOCATE(aqsmtemp(ncouls,number_bands))
      aqsmtemp = (0.5D0,0.5D0)
      aqsntemp = (0.5D0,0.5D0)

      write(6,*)"allocating I_eps*_array, size = ",((2.0*ncouls)*ngpown*nFreq*16)/1024**2," Mbytes"
      ALLOCATE(I_epsR_array(ncouls,ngpown,nFreq))
      I_epsR_array = (0.5D0,0.5D0)
      ALLOCATE(I_epsA_array(ncouls,ngpown,nFreq))
      I_epsA_array = (0.5D0,-0.5D0)

      write(6,*)"allocating matn*, size = ",2*ncouls*ngpown*16/1024," kbytes"
      ALLOCATE(matngmatmgpD(ncouls,ngpown))
      ALLOCATE(matngpmatmgD(ncouls,ngpown))

      ALLOCATE(inv_igp_index(ngpown))
      ALLOCATE(indinv(ncouls))

      ALLOCATE(dFreqGrid(nFreq))
      dw = 0D0
      do ijk = 1, nFreq
        dFreqGrid(ijk) = dw
        dw = dw + 2D0
      enddo

      ALLOCATE(pref(nFreq))
      do ifreq=1,nFreq
        if (ifreq .lt. nFreq) then
          pref(ifreq)=(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))/3.14d0
        else
          pref(ifreq)=pref(ifreq-1)
        endif
      enddo
      pref(1)=pref(1)*0.5d0
      pref(nFreq)=pref(nFreq)*0.5d0

      ALLOCATE(dFreqBrd(nFreq))
      dFreqBrd = (0D0,0.1D0)

      ALLOCATE(asxDtemp(nfreqeval))
      asxDtemp = 0D0
      ALLOCATE(achDtemp(nfreqeval))
      achDtemp = 0D0
      ALLOCATE(achDtemp_cor(nfreqeval))
      achDtemp_cor = 0D0
      ALLOCATE(achDtemp_corb(nfreqeval))
      achDtemp_corb = 0D0
      ALLOCATE(ach2Dtemp(nfreqeval))
      ach2Dtemp = 0D0
      ALLOCATE(schDi(nfreqeval))
      schDi=0D0
      ALLOCATE(schDi_cor(nfreqeval))
      schDi_cor=0D0
      ALLOCATE(schDi_corb(nfreqeval))
      schDi_corb=0D0
      ALLOCATE(sch2Di(nfreqeval))
      sch2Di=0D0
      ALLOCATE(ssxDi(nfreqeval))
      ssxDi=0D0

      ALLOCATE(wxi(nfreqeval))
      wxi=0D0
      ALLOCATE(schDt_matrix(nFreq,number_bands))
      schDt_matrix = 0.

      do ig = 1, ngpown
        inv_igp_index(ig) = ig
      enddo

! Try this identity and random
      do ig = 1, ncouls
        indinv(ig) = ig
      enddo

      tol = 1D-6
      limitone=1D0/(tol*4D0)
      limittwo=0.5d0**2

      write(6,*) "Starting loop"

      call timget(starttime)

      do n1=1,number_bands

! n1true = "True" band index of the band n1 w.r.t. all bands

        !!! n1true = peinf%indext_dist(n1,ipe)

! energy of the |n1,k-q> state

        e_n1kq = ekq(n1,1)

! occupation of the |n1,k-q> state

        flag_occ = (n1.le.nvband)

        !!!tempval=abs(e_n1kq-efermi)
        !!!if (tempval .lt. tol) then
        !!!  occ = 0.5d0 ! Fermi-Dirac distribution = 1/2 at Fermi level
        !!!else
          occ = 1.0d0
        !!!endif

! JRD: compute the static CH for the static remainder

        iwlda = 1

!         do iw=1,nfreqeval
!           wx = freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep
!           wxi(iw) = wx
!         enddo

! JRD compute the static CH for the static remainder

        call timget(starttime_stat)

! !$OMP PARALLEL do private (my_igp,igp,indigp,igmax,mygpvar1,ig)
!         do my_igp = 1, ngpown
!           indigp = inv_igp_index(my_igp)
!           igp = indinv(indigp)
! 
!           if (igp .gt. ncouls .or. igp .le. 0) cycle
! 
!           igmax=ncouls
! 
!           mygpvar1 = CONJG(aqsmtemp(igp,n1))
! 
!           do ig = 1, igmax
!             matngmatmgpD(ig,my_igp) = aqsntemp(ig,n1) * CONJG(aqsmtemp(igp,n1))
! !             matngmatmgpD(ig,my_igp) = aqsntemp(ig,n1) * mygpvar1
!           enddo
!         enddo
! !$OMP END PARALLEL DO

        !if (exact_ch.eq.1) then
!$OMP PARALLEL do private (my_igp,igp,indigp,igmax,ig, &
!$OMP                      schsDtemp ) reduction(+:achsDtemp)
          do my_igp = 1, ngpown
            indigp = inv_igp_index(my_igp)
            igp = indinv(indigp)

            if (igp .gt. ncouls .or. igp .le. 0) cycle

            igmax=ncouls

            schsDtemp = 0D0
            do ig = 1, igmax
              schsDtemp = schsDtemp-aqsntemp(ig,n1) * CONJG(aqsmtemp(igp,n1))*I_epsR_array(ig,my_igp,1)
            enddo
            achsDtemp = achsDtemp + schsDtemp*vcoul(igp)*0.5D0
          enddo
!$OMP END PARALLEL DO
        !endif

        call timget(endtime_stat)
        time_stat = time_stat + endtime_stat - starttime_stat

        ssxDi = (0D0,0D0)

        call timget(starttime_sx)

! JRD: Don`t want to thread here, nfreqeval could be small
        do iw=1,nfreqeval
            
          wx = freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep

! SX and CH terms: equation (1.42) of Catalin`s thesis
! Note the negative sign in I_epsRggp and I_epsAggp

          if (flag_occ) then

            if(wx.ge.0.0d0) then
              ifreq=0
              do ijk = 1, nFreq-1
                if (wx .ge. dFreqGrid(ijk) .and. wx .lt. dFreqGrid(ijk+1)) then
                  ifreq=ijk
                endif
              enddo
              if (ifreq .eq. 0) then
                ifreq = nFreq+3 ! three is for luck
              endif
            else
              ifreq=0
              do ijk = 1, nFreq-1
                if (-wx .ge. dFreqGrid(ijk) .and. -wx .lt. dFreqGrid(ijk+1)) then
                  ifreq=ijk
                endif
              enddo
              if (ifreq .eq. 0) then
                ifreq = nFreq+3 ! three is for luck
              endif
            endif
              
            if(ifreq.ge.nFreq) then
              !write(0,777) n1,e_lk,e_n1kq,wx
              ifreq=nFreq-1
            endif
777         format(1x,"WARNING: The real frequency range is too small." &
              ,/,3x,"n1 =",i5,1x,"E_l =",f8.3,1x,"E_n1" &
              ,1x,"=",f8.3,1x,"wx =",f8.3,1x)
              
            if(wx.ge.0.d0) then

              fact1 = (dFreqGrid(ifreq+1)-wx)/(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))
              fact2 = (wx-dFreqGrid(ifreq))/(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))

              ssxDittt = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,igmax,ssxDitt,ig, &
!$OMP                       ssxDit) reduction(+:ssxDittt)
              do my_igp = 1, ngpown
                indigp = inv_igp_index(my_igp)
                igp = indinv(indigp)

                if (igp .gt. ncouls .or. igp .le. 0) cycle

                igmax=ncouls

                ssxDitt = (0D0,0D0)
                do ig = 1, igmax
                  ssxDit=I_epsR_array(ig,my_igp,ifreq)*fact1 + &
                  I_epsR_array(ig,my_igp,ifreq+1)*fact2 
                  !ssxDit=I_epsR_array(ifreq,ig,my_igp)*fact1 + &
                  !I_epsR_array(ifreq+1,ig,my_igp)*fact2 
 
                  ssxDitt = ssxDitt + aqsntemp(ig,n1) * CONJG(aqsmtemp(igp,n1))*ssxDit
                enddo
                ssxDittt = ssxDittt + ssxDitt*vcoul(igp)
              enddo
!$OMP END PARALLEL DO

              ssxDi(iw) = ssxDi(iw) + ssxDittt

            else

              fact1 = (dFreqGrid(ifreq+1)+wx)/(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))
              fact2 = (-dFreqGrid(ifreq)-wx)/(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))

              ssxDittt = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,igmax,ssxDitt,ig, &
!$OMP                       ssxDit) reduction(+:ssxDittt)
              do my_igp = 1, ngpown
                indigp = inv_igp_index(my_igp)
                igp = indinv(indigp)

                if (igp .gt. ncouls .or. igp .le. 0) cycle

                igmax=ncouls

                ssxDitt = (0D0,0D0)
                do ig = 1, igmax
                  ssxDit=I_epsA_array(ig,my_igp,ifreq)*fact1+ &
                    I_epsA_array(ig,my_igp,ifreq+1)*fact2
                  !ssxDit=I_epsA_array(ifreq,ig,my_igp)*fact1+ &
                  !  I_epsA_array(ifreq+1,ig,my_igp)*fact2

                  ssxDitt = ssxDitt + aqsntemp(ig,n1) * CONJG(aqsmtemp(igp,n1))*ssxDit
                enddo
                ssxDittt = ssxDittt + ssxDitt*vcoul(igp)
              enddo                
!$OMP END PARALLEL DO

              ssxDi(iw) = ssxDi(iw) + ssxDittt

            endif
          endif
          if (flag_occ) then
            asxDtemp(iw) = asxDtemp(iw) + ssxDi(iw)*occ
          endif
        enddo

        call timget(endtime_sx)
        time_sx = time_sx + endtime_sx - starttime_sx
      enddo

! JRD: Now do CH term
          
!         ALLOCATE(schDt_array(nFreq),schDt_matrix(nFreq,))

        call timget(starttime_ch)
 
!         schdt_array = 0D0
!$OMP PARALLEL do private (my_igp,igp,indigp,igmax,ig,schDtt,I_epsRggp_int, &
!$OMP                      I_epsAggp_int,schD,schDt,ifreq)

! The following omp directive is in error.  There is not a reduction on schdt_array with
! respect to the omp loop, in this case the ifreq loop.  Each omp iteration has a unique
! location so it is not a reduction.
! ! ! !$OMP                      I_epsAggp_int,schD,schDt,ifreq) reduction(+:schdt_array)
        do ifreq=1,nFreq

!             schDt = (0D0,0D0)

            do my_igp = 1, ngpown
              indigp = inv_igp_index(my_igp)
              igp = indinv(indigp)

              if (igp .gt. ncouls .or. igp .le. 0) cycle

              igmax=ncouls

! JRD: The below loop is performance critical

      do n1=1,number_bands
              schDtt = (0D0,0D0)
              do ig = 1, igmax

                I_epsRggp_int = I_epsR_array(ig,my_igp,ifreq)

                I_epsAggp_int = I_epsA_array(ig,my_igp,ifreq)

                ! for G,G` components
                schD=I_epsRggp_int-I_epsAggp_int

                ! for G`,G components
                schDtt = schDtt +aqsntemp(ig,n1) * CONJG(aqsmtemp(igp,n1))*schD
              enddo
              !schDt = schDt + schDtt * vcoul(igp)
              schdt_matrix(ifreq,n1) = schdt_matrix(ifreq,n1) + schDtt
      enddo
            enddo

! XXX: Threads could be stomping on each-other`s cache over this... try reduction?
!            schdt_array(ifreq) = schdt_array(ifreq) + schDt

        enddo
!$OMP END PARALLEL DO

        call timget(endtime_ch)
        time_cha = time_cha + endtime_ch - starttime_ch

      do n1=1,number_bands
        flag_occ = (n1.le.nvband)
        occ = 1.0d0
        call timget(starttime_ch)
        schDi = (0D0,0D0)
        schDi_cor = (0D0,0D0)
        schDi_corb = (0D0,0D0)
        sch2Di = (0D0,0D0)

!$OMP PARALLEL do private (ifreq,schDt,cedifft_zb,cedifft_coh,cedifft_cor, &
!$OMP                      cedifft_zb_right,cedifft_zb_left,schDt_right,schDt_left, &
!$OMP                      schDt_avg,schDt_lin,schDt_lin2,intfact,iw, &
!$OMP                      schDt_lin3) reduction(+:schDi,schDi_corb,schDi_cor,sch2Di)
        do ifreq=1,nFreq

            schDt = schdt_matrix(ifreq,n1)

            cedifft_zb = dFreqGrid(ifreq)
            cedifft_coh = CMPLX(cedifft_zb,0D0)- dFreqBrd(ifreq)

            if( flag_occ) then 
              cedifft_cor = -1.0d0*CMPLX(cedifft_zb,0D0) - dFreqBrd(ifreq)
            else
              cedifft_cor = CMPLX(cedifft_zb,0D0) - dFreqBrd(ifreq)
            endif

            if (ifreq .ne. 1) then 
              cedifft_zb_right = cedifft_zb
              cedifft_zb_left = dFreqGrid(ifreq-1)
              schDt_right = schDt
              schDt_left = schdt_matrix(ifreq-1,n1)
              schDt_avg = 0.5D0 * ( schDt_right + schDt_left )
              schDt_lin = schDt_right - schDt_left
              schDt_lin2 = schDt_lin/(cedifft_zb_right-cedifft_zb_left)
            endif

! The below two lines are for sigma1 and sigma3
            if (ifreq .ne. nFreq) then
              do iw = 1, nfreqeval
                wx = freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep
                schDi(iw) = schDi(iw) - CMPLX(0.d0,pref(ifreq)) * schDt / ( wx-cedifft_coh)
                schDi_corb(iw) = schDi_corb(iw) - CMPLX(0.d0,pref(ifreq)) * schDt / ( wx-cedifft_cor)
              enddo
            endif
            if (ifreq .ne. 1) then
              do iw = 1, nfreqeval
!These lines are for sigma2
                intfact=abs((freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep-cedifft_zb_right)/(freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep-cedifft_zb_left))
                if (intfact .lt. 1d-4) intfact = 1d-4
                if (intfact .gt. 1d4) intfact = 1d4
                intfact = -log(intfact)
                sch2Di(iw) = sch2Di(iw) - CMPLX(0.d0,pref_zb) * schDt_avg * intfact
!These lines are for sigma4
                if (flag_occ) then
                  intfact=abs((freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep+cedifft_zb_right)/(freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep+cedifft_zb_left))
                  if (intfact .lt. 1d-4) intfact = 1d-4
                  if (intfact .gt. 1d4) intfact = 1d4
                  intfact = log(intfact)
                  schDt_lin3 = (schDt_left + schDt_lin2*(-freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep-cedifft_zb_left))*intfact
                else 
                  schDt_lin3 = (schDt_left + schDt_lin2*(freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep-cedifft_zb_left))*intfact
                endif
                schDt_lin3 = schDt_lin3 + schDt_lin
                schDi_cor(iw) = schDi_cor(iw) - CMPLX(0.d0,pref_zb) * schDt_lin3
              enddo
            endif
        enddo
!$OMP END PARALLEL DO

!         DEALLOCATE(schDt_array)

! JRD: Compute Sigma2 and Sigma4 delta function contributions

        call timget(endtime_ch)
        time_chb = time_chb + endtime_ch - starttime_ch
        call timget(starttime_ch)

! JRD: This can be often (but not always) small, so we don't thread over it. Maybe smarter to dynamically whether to thread this loop
! or below loop, or just use OpenMP nested loop support
        do iw = 1, nfreqeval
          wx = freqevalmin - ekq(n1,1) + (iw-1)*freqevalstep
          if(wx .ge. 0.0d0) then
            ifreq=0
            do ijk = 1, nFreq-1
              if (wx .ge. dFreqGrid(ijk) .and. wx .lt. dFreqGrid(ijk+1)) then
                ifreq=ijk
              endif
            enddo
            if (ifreq .eq. 0) then
              ifreq=nFreq-1
            endif

            fact1=-0.5D0*(dFreqGrid(ifreq+1)-wx)/(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))
            fact2=-0.5D0*(wx-dFreqGrid(ifreq))/(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))

            schDttt = 0D0
            schDttt_cor = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,igmax,ig, &
!$OMP                      sch2Dt,sch2Dtt) reduction(+:schDttt,schDttt_cor) 
            do my_igp = 1, ngpown
              indigp = inv_igp_index(my_igp)
              igp = indinv(indigp)

              if (igp .gt. ncouls .or. igp .le. 0) cycle

              igmax=ncouls

              sch2Dtt = (0D0,0D0)
              do ig = 1, igmax
                sch2Dt=(I_epsR_array(ig,my_igp,ifreq)-I_epsA_array(ig,my_igp,ifreq))*fact1 + &
                       (I_epsR_array(ig,my_igp,ifreq+1)-I_epsA_array(ig,my_igp,ifreq+1))*fact2
                sch2Dtt = sch2Dtt + aqsntemp(ig,n1) * CONJG(aqsmtemp(igp,n1))*sch2Dt
              enddo
              schDttt = schDttt + sch2Dtt*vcoul(igp)
              if (flag_occ) then
              else
                schDttt_cor = schDttt_cor + sch2Dtt*vcoul(igp)
              endif
            enddo
!$OMP END PARALLEL DO

            sch2Di(iw) = sch2Di(iw) + schDttt
            schDi_cor(iw) = schDi_cor(iw) + schDttt_cor
          else if (flag_occ) then
            wx=-wx
            ifreq=0
            do ijk = 1, nFreq-1
              if (wx .ge. dFreqGrid(ijk) .and. wx .lt. dFreqGrid(ijk+1)) then
                ifreq=ijk
              endif
            enddo
            if (ifreq .eq. 0) then
              ifreq=nFreq-1
            endif

            fact1=(dFreqGrid(ifreq+1)-wx)/(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))
            fact2=(wx-dFreqGrid(ifreq))/(dFreqGrid(ifreq+1)-dFreqGrid(ifreq))

            schDttt_cor = 0D0

!$OMP PARALLEL do private (my_igp,igp,indigp,igmax,ig, &
!$OMP                      sch2Dt,sch2Dtt) reduction(+:schDttt_cor) 
            do my_igp = 1, ngpown
              indigp = inv_igp_index(my_igp)
              igp = indinv(indigp)

              if (igp .gt. ncouls .or. igp .le. 0) cycle

              igmax=ncouls

              sch2Dtt = (0D0,0D0)
              do ig = 1, igmax
                sch2Dt=-0.5D0*((I_epsR_array(ig,my_igp,ifreq)-I_epsA_array(ig,my_igp,ifreq))*fact1 + &
                       (I_epsR_array(ig,my_igp,ifreq+1)-I_epsA_array(ig,my_igp,ifreq+1))*fact2)
                sch2Dtt = sch2Dtt + aqsntemp(ig,n1) * CONJG(aqsmtemp(igp,n1))*sch2Dt
              enddo
              schDttt_cor = schDttt_cor + sch2Dtt*vcoul(igp)
            enddo
!$OMP END PARALLEL DO
            schDi_cor(iw) = schDi_cor(iw) + schDttt_cor
          endif
        enddo

        call timget(endtime_ch)
        time_chc = time_chc + endtime_ch - starttime_ch

        do iw = 1, nfreqeval
            
!           if (flag_occ) then
!             asxDtemp(iw) = asxDtemp(iw) + ssxDi(iw)*occ
!           endif
            
          achDtemp(iw) = achDtemp(iw) + schDi(iw)
          achDtemp_cor(iw) = achDtemp_cor(iw) + schDi_cor(iw)
          achDtemp_corb(iw) = achDtemp_corb(iw) + schDi_corb(iw)
          ach2Dtemp(iw) = ach2Dtemp(iw) + sch2Di(iw)

        enddo ! over iw

      enddo ! over ipe bands (n1)

      call timget(endtime)

      DEALLOCATE(vcoul)
      DEALLOCATE(aqsntemp)
      DEALLOCATE(aqsmtemp)
      DEALLOCATE(inv_igp_index)
      DEALLOCATE(indinv)
      DEALLOCATE(I_epsR_array)
      DEALLOCATE(I_epsA_array)
      DEALLOCATE(ekq)
      DEALLOCATE(matngmatmgpD)
      DEALLOCATE(matngpmatmgD)
      DEALLOCATE(achDtemp_corb)
      DEALLOCATE(ach2Dtemp)
      DEALLOCATE(schDi)
      DEALLOCATE(schDi_cor)
      DEALLOCATE(schDi_corb)
      DEALLOCATE(sch2Di)
      DEALLOCATE(ssxDi)
      DEALLOCATE(wxi)
      DEALLOCATE(dFreqBrd)
      DEALLOCATE(dFreqGrid)
      DEALLOCATE(schDt_matrix)
      

      write(6,*) "Runtime:", endtime-starttime
      write(6,*) "Runtime STAT:", time_stat
      write(6,*) "Runtime SX:", time_sx
      write(6,*) "Runtime CH:", time_cha, time_chb, time_chc
      write(6,*) "Answer:", achDtemp_cor(1), achDtemp(1), asxDtemp(1)

      DEALLOCATE(asxDtemp)
      DEALLOCATE(achDtemp)
      DEALLOCATE(achDtemp_cor)

end program

subroutine timget(wall)
  real(kind(1.0d0)) :: wall

  integer :: values(8)

  call date_and_time(VALUES=values)
  wall=((values(3)*24.0d0+values(5))*60.0d0 &
    +values(6))*60.0d0+values(7)+values(8)*1.0d-3

  return
end subroutine timget
