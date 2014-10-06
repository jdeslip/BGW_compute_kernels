!=============================================================================================
!
! This Kernel Represents the GPP Self-Energy Summations in BerkeleyGW
!
! Run like:
!
!  gppkernel.x <number_bands> <number_valence_bands> <number_plane_waves> <number_plane_waves_per_mpi_task>
!
! Example run:
!
!  gppkernel.x 8 2 10000 600
!
!==============================================================================================


program gppkernel

      implicit none

      integer :: NTHREADS,TID,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
      integer :: n1
      integer :: ispin, iw
      integer :: number_bands,nvband,ncouls,ngpown
      integer :: nstart, nend
      integer :: my_igp, indigp, ig, igp, igmax
      integer :: ggpsum
      integer, allocatable :: indinv(:), inv_igp_index(:)
      complex(kind((1.0d0,1.0d0))) :: achstemp,achxtemp,matngmatmgp,matngpmatmg,wdiff,mygpvar1,mygpvar2,schstemp,schs,wx_array(3),sch,ssx,ssxt,scht
      complex(kind((1.0d0,1.0d0))) :: wtilde,halfinvwtilde,wtilde2,Omega2
      complex(kind((1.0d0,1.0d0))), allocatable :: aqsmtemp(:,:), aqsntemp(:,:), I_eps_array(:,:), acht_n1_loc(:)
      complex(kind((1.0d0,1.0d0))), allocatable :: asxtemp(:),achtemp(:),ssx_array(:),sch_array(:),ssxa(:),scha(:),wtilde_array(:,:)
      complex(kind((1.0d0,1.0d0))) :: delw, cden
      real(kind(1.0d0)) :: e_n1kq, e_lk, dw, delwr, wdiffr
      real(kind(1.0d0)) :: tol, gamma, sexcut
      real(kind(1.0d0)) :: wx,occ,occfact,tempval
      real(kind(1.0d0)) :: wxt,rden,ssxcutoff,delw2
      real(kind(1.0d0)), allocatable :: ekq(:,:), vcoul(:)
      real(kind(1.0d0)) :: limitone,limittwo
      real(kind(1.0d0)) :: starttime, endtime
      real(kind(1.0d0)) :: starttime_stat, endtime_stat
      real(kind(1.0d0)) :: time_stat
      real(kind(1.0d0)) :: starttime_dyn, endtime_dyn
      real(kind(1.0d0)) :: time_dyn
      logical :: flag_occ
      CHARACTER(len=32) :: arg

      time_stat = 0D0
      time_dyn = 0D0

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

! ngpown = ncouls / (number of mpi tasks)

      e_lk = 10D0
      dw = 1D0
      nstart = 1
      nend = 3
      ggpsum = 2

      ALLOCATE(vcoul(ncouls))
      vcoul = 1D0

      ALLOCATE(asxtemp(nend-nstart+1))
      asxtemp = 0D0

      ALLOCATE(achtemp(nend-nstart+1))
      achtemp = 0D0

      ALLOCATE(ekq(number_bands,1))
      ekq = 6D0 

      ALLOCATE(wtilde_array(ncouls,ngpown))
      wtilde_array = (0.5d0,0.5d0)

      ALLOCATE(aqsntemp(ncouls,number_bands))
      ALLOCATE(aqsmtemp(ncouls,number_bands))
      aqsmtemp = (0.5D0,0.5D0)
      aqsntemp = (0.5D0,0.5D0)

      ALLOCATE(I_eps_array(ncouls,ngpown))
      I_eps_array = (0.5D0,0.5D0)

      ALLOCATE(inv_igp_index(ngpown))
      ALLOCATE(indinv(ncouls))

      ALLOCATE(acht_n1_loc(number_bands))

      do ig = 1, ngpown
        inv_igp_index(ig) = ig
      enddo

! Try this as identity and random
! Identity for now
      do ig = 1, ncouls
        indinv(ig) = ig
      enddo

      tol = 1D-6
      gamma = 0.5d0
      sexcut=4.0d0
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

        !!!tempval=abs(e_n1kq-sig%efermi)
        !!!if (tempval .lt. tol) then
        !!!  occ = 0.5d0 ! Fermi-Dirac distribution = 1/2 at Fermi level
        !!!else
          occ = 1.0d0
        !!!endif

! JRD: compute the static CH for the static remainder

        call timget(starttime_stat)
 
        !!!if (sig%exact_ch.eq.1) then
!$OMP PARALLEL do private (my_igp,igp,indigp,igmax,mygpvar1,mygpvar2,ig,schs, &
!$OMP                       matngmatmgp,matngpmatmg,schstemp) reduction(+:achstemp)
          do my_igp = 1, ngpown

            ! This index could be bad

            indigp = inv_igp_index(my_igp)
            igp = indinv(indigp)

            if (igp .gt. ncouls .or. igp .le. 0) cycle

            if (ggpsum.eq.1) then
              igmax=igp
            else
              igmax=ncouls
            endif

            mygpvar1 = CONJG(aqsmtemp(igp,n1))
            mygpvar2 = aqsntemp(igp,n1)

            schstemp = 0D0

! We do two loops here for performance. Don`t want to evaluate if statements inside loop
! at every iteration

            if (ggpsum.eq.1) then
              do ig = 1, igmax - 1
                schs=-I_eps_array(ig,my_igp)
! JRD: Cycle bad for vectorization.
! I_eps_array is already set to zero above for these ig,igp
!                if (abs(schs).lt.tol) cycle
                matngmatmgp = aqsntemp(ig,n1) * mygpvar1
                matngpmatmg = CONJG(aqsmtemp(ig,n1)) * mygpvar2
                schstemp = schstemp + matngmatmgp*schs + matngpmatmg*CONJG(schs)
              enddo
              ig = igp
              schs=-I_eps_array(ig,my_igp)
              matngmatmgp = aqsntemp(ig,n1) * mygpvar1
              if (abs(schs).gt.tol) schstemp = schstemp + matngmatmgp*schs
            else
              do ig = 1, igmax
                !schs=-I_eps_array(ig,my_igp)
! JRD: Cycle bad for vectorization.
! I_eps_array is already set to zero above for these ig,igp
!             if (abs(schs).lt.tol) cycle
                !matngmatmgp = aqsntemp(ig,n1) * mygpvar1
                schstemp = schstemp - aqsntemp(ig,n1) * I_eps_array(ig,my_igp) * mygpvar1
                !schstemp = schstemp + matngmatmgp * schs
              enddo
            endif

            achstemp = achstemp + schstemp*vcoul(igp)*0.5d0
          enddo
!$OMP END PARALLEL DO
        !!!endif ! sig%exact_ch.eq.1

        call timget(endtime_stat)

        time_stat = time_stat + endtime_stat - starttime_stat

        do iw=nstart,nend
          wx_array(iw) = e_lk - e_n1kq + dw*(iw-2)
          if (abs(wx_array(iw)) .lt. tol) wx_array(iw) = tol
        enddo

        call timget(starttime_dyn)

! JRD: This Loop is Performance critical. Make Sure you don`t mess it up

!$OMP PARALLEL private (my_igp,igp,indigp,igmax,mygpvar1,mygpvar2,ssx_array,sch_array,ig, &
!$OMP                      wtilde,wtilde2,halfinvwtilde,ssxcutoff,matngmatmgp,matngpmatmg,sch,ssx, &
!$OMP                      iw,delw,delw2,Omega2,scht,ssxt,wxt, &
!$OMP                      rden,cden,ssxa,scha,delwr,wdiffr,occfact)

        ALLOCATE(ssx_array(3))
        ALLOCATE(sch_array(3))

        ALLOCATE(ssxa(ncouls))
        ALLOCATE(scha(ncouls))

!$OMP DO reduction(+:asxtemp,acht_n1_loc,achtemp)
        do my_igp = 1, ngpown

          indigp = inv_igp_index(my_igp)
          igp = indinv(indigp)

          if (igp .gt. ncouls .or. igp .le. 0) cycle

          if (ggpsum.eq.1) then
            igmax=igp
          else
            igmax=ncouls
          endif

          ssx_array = 0D0
          sch_array = 0D0

          mygpvar1 = CONJG(aqsmtemp(igp,n1))
          mygpvar2 = aqsntemp(igp,n1)

          if (flag_occ) then

            do iw=nstart,nend

              scht=0D0
              ssxt=0D0
              wxt = wx_array(iw)

              if (ggpsum.eq.1) then

                do ig = 1, igmax
                  wtilde = wtilde_array(ig,my_igp)
                  wtilde2 = wtilde**2
                  Omega2 = wtilde2 * I_eps_array(ig,my_igp)

! Cycle bad for vectorization. ggpsum=1 slow anyway
                  if (abs(Omega2) .lt. tol) cycle

                  matngmatmgp = aqsntemp(ig,n1) * mygpvar1
! JRD: This breaks vectorization but ggpsum=1 is slow for other reasons already
                  if (ig .ne. igp) matngpmatmg = CONJG(aqsmtemp(ig,n1)) * mygpvar2

                  halfinvwtilde = 0.5d0/wtilde
                  delw = (wxt - wtilde) * halfinvwtilde
                  delw2 = abs(delw)**2
                  if (abs(wxt-wtilde).lt.gamma .or. delw2.lt.tol) then
                    sch = 0.0d0
                    if (abs(wtilde) .gt. tol) then
                      ssx = -Omega2 / (4.0d0 * wtilde2 * (1.0d0 + delw))
                    else
                      ssx = 0D0
                    endif
                  else
                    sch = wtilde * I_eps_array(ig,my_igp) / (wxt - wtilde)
                    ssx = Omega2 / (wxt**2 - wtilde2)
                  endif

! JRD: Bad for vectorization
                  ssxcutoff = sexcut*abs(I_eps_array(ig,my_igp))
                  if (abs(ssx) .gt. ssxcutoff .and. wxt .lt. 0.0d0) ssx=0.0d0

                  if (ig .ne. igp) then
                    ssxa(ig) = matngmatmgp*ssx + matngpmatmg*CONJG(ssx)
                    scha(ig) = matngmatmgp*sch + matngpmatmg*CONJG(sch)
                  else
                    ssxa(ig) = matngmatmgp*ssx
                    scha(ig) = matngmatmgp*sch
                  endif

                  ssxt = ssxt + ssxa(ig)
                  scht = scht + scha(ig)

                enddo ! loop over g

              else

                do ig = 1, igmax

                  wtilde = wtilde_array(ig,my_igp)
                  wtilde2 = wtilde**2
                  Omega2 = wtilde2 * I_eps_array(ig,my_igp)

! Cycle bad for vectorization. Not needed wtilde is zero
                  !if (abs(Omega2) .lt. tol) cycle

                  matngmatmgp = aqsntemp(ig,n1) * mygpvar1

                  wdiff = wxt - wtilde

                  cden = wdiff
                  rden = cden * CONJG(cden)
                  rden = 1D0 / rden
                  delw = wtilde * CONJG(cden) * rden
                  delwr = delw*CONJG(delw)
                  wdiffr = wdiff*CONJG(wdiff)

! This Practice is bad for vectorization and understanding of the output.
! JRD: Complex division is hard to vectorize. So, we help the compiler.
                  if (wdiffr.gt.limittwo .and. delwr.lt.limitone) then
                    sch = delw * I_eps_array(ig,my_igp)
                    cden = wxt**2 - wtilde2
                    rden = cden*CONJG(cden)
                    rden = 1D0 / rden
                    ssx = Omega2 * CONJG(cden) * rden
                  else if ( delwr .gt. tol) then
                    sch = 0.0d0
                    cden = (4.0d0 * wtilde2 * (delw + 0.5D0 ))
                    rden = cden*CONJG(cden)
                    rden = 1D0 / rden
                    ssx = -Omega2 * CONJG(cden) * rden * delw
                  else
                    sch = 0.0d0
                    ssx = 0.0d0
                  endif

! JRD: Breaks vectorization. But, I will have to fix later because
! leaving it out breaks GSM example.
                  ssxcutoff = sexcut*abs(I_eps_array(ig,my_igp))
                  if (abs(ssx) .gt. ssxcutoff .and. wxt .lt. 0.0d0) ssx=0.0d0

                  ssxa(ig) = matngmatmgp*ssx
                  scha(ig) = matngmatmgp*sch

                  ssxt = ssxt + ssxa(ig)
                  scht = scht + scha(ig)

                enddo ! loop over g

              endif

              ssx_array(iw) = ssx_array(iw) + ssxt
              sch_array(iw) = sch_array(iw) + 0.5D0*scht

            enddo

          else

            do iw=nstart,nend

              scht=0D0
              ssxt=0D0
              wxt = wx_array(iw)

              if (ggpsum.eq.1) then

                do ig = 1, igmax

                  wtilde = wtilde_array(ig,my_igp)

! Cycle bad for vectorization. ggpsum=1 slow anyway
                  if (abs((wtilde**2) * I_eps_array(ig,my_igp)) .lt. tol) cycle

                  matngmatmgp = aqsntemp(ig,n1) * mygpvar1
! JRD: If statement breaks vectorization but gppsum=1 is slow anyway
                  if (ig .ne. igp) matngpmatmg = CONJG(aqsmtemp(ig,n1)) * mygpvar2

                  halfinvwtilde = 0.5d0/wtilde
                  delw = (wx_array(iw) - wtilde) * halfinvwtilde
                  delw2 = abs(delw)**2

                  if (abs(wxt-wtilde).lt.gamma .or. delw2.lt.tol) then
                    sch = 0.0d0
                  else
                    sch = wtilde_array(ig,my_igp) * I_eps_array(ig,my_igp) / (wxt - wtilde_array(ig,my_igp))
                  endif

                  if (ig .ne. igp) then
                    scha(ig) = matngmatmgp*sch + matngpmatmg*CONJG(sch)
                  else
                    scha(ig) = matngmatmgp*sch
                  endif

                  scht = scht + scha(ig)

                enddo ! loop over g

              else

                do ig = 1, igmax

                  wdiff = wxt - wtilde_array(ig,my_igp)

                  cden = wdiff
                  rden = cden * CONJG(cden)
                  rden = 1D0 / rden
                  delw = wtilde_array(ig,my_igp) * CONJG(cden) * rden
                  delwr = delw*CONJG(delw)
                  wdiffr = wdiff*CONJG(wdiff)

! JRD: Complex division is hard to vectorize. So, we help the compiler.
                  scha(ig) = mygpvar1 * aqsntemp(ig,n1) * delw * I_eps_array(ig,my_igp)

! JRD: This if is OK for vectorization
                  if (wdiffr.gt.limittwo .and. delwr.lt.limitone) then
                    scht = scht + scha(ig)
                  endif

                enddo ! loop over g

              endif

              sch_array(iw) = sch_array(iw) + 0.5D0*scht

!-----------------------
! JRD: Compute GPP Error...
! GPP Model Error Estimate

            enddo

          endif

! If a valence band, then accumulate SX contribution.

          if (flag_occ) then
            do iw=nstart,nend
              asxtemp(iw) = asxtemp(iw) - ssx_array(iw) * occ * vcoul(igp)
            enddo
          endif

          do iw=nstart,nend
            achtemp(iw) = achtemp(iw) + sch_array(iw) * vcoul(igp)
          enddo

! Logging CH convergence.

          acht_n1_loc(n1) = acht_n1_loc(n1) + sch_array(2) * vcoul(igp)

        enddo ! igp
!$OMP END DO

        DEALLOCATE(ssx_array)
        DEALLOCATE(sch_array)
        DEALLOCATE(ssxa)
        DEALLOCATE(scha)

!$OMP END PARALLEL
        call timget(endtime_dyn)
        time_dyn = time_dyn + endtime_dyn - starttime_dyn

      enddo ! over ipe bands (n1)

      call timget(endtime)

      DEALLOCATE(vcoul)
      DEALLOCATE(aqsntemp)
      DEALLOCATE(aqsmtemp)
      DEALLOCATE(inv_igp_index)
      DEALLOCATE(indinv)
      DEALLOCATE(I_eps_array)
      DEALLOCATE(acht_n1_loc)
      DEALLOCATE(wtilde_array)
      DEALLOCATE(ekq)

      write(6,*) "Runtime:", endtime-starttime
      write(6,*) "Runtime Stat:", time_stat
      write(6,*) "Runtime Dyn:", time_dyn
      write(6,*) "Answer:",achtemp(2)

      DEALLOCATE(achtemp)
      DEALLOCATE(asxtemp)

end program

subroutine timget(wall)
  real(kind(1.0d0)) :: wall

  integer :: values(8)

  call date_and_time(VALUES=values)
  wall=((values(3)*24.0d0+values(5))*60.0d0 &
    +values(6))*60.0d0+values(7)+values(8)*1.0d-3

  return
end subroutine timget
