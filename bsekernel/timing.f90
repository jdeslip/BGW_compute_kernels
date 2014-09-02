module timing_m

  use common_m

  implicit none

  private

  public ::    &
    date_time, &
    timget,    &
    timacc

  !> mtim determines the maximum number of "timing slots" available
  integer, parameter, private :: mtim=100
  real(DP), private, save :: acctim(2,mtim),tzero(2,mtim)
  integer, private, save :: ncount(mtim)

contains

  subroutine date_time(bdate,btime)
    character, intent(out) :: bdate*11,btime*14
    
    integer :: lmonth
    integer :: idate (8)
    character :: day*2,year*4
    character :: adate*8,atime*10,azone*5
    character :: hour*2,min*2,sec*2
    character*3 :: month(12)
      
    DATA month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', &
      'Oct','Nov','Dec'/

    call date_and_time(adate,atime,azone,idate)
    
    read(adate,101) year,lmonth,day
101 format(a4,i2,a2)
    write(bdate,102) day,month(lmonth),year
102 format(a2,'-',a3,'-',a4)
    read(atime,201) hour,min,sec
201 format(a2,a2,a2,a4)
    write(btime,202) hour,min,sec,azone
202 format(a2,':',a2,':',a2,1x,a5)

    return    
  end subroutine date_time

!================================================================================

  subroutine timget(cpu,wall)
    real(DP), intent(out) :: cpu,wall
    
    integer :: values(8)

    ! no push_sub, called too frequently

!XXX
!intel
!    cpu=mclock()*1.0d-3
!cray
    call cpu_time(cpu)

    call date_and_time(VALUES=values)
    wall=((values(3)*24.0d0+values(5))*60.0d0 &
      +values(6))*60.0d0+values(7)+values(8)*1.0d-3

    return
  end subroutine timget

!================================================================================

  subroutine timacc(n,option,tottim,nslices)
    integer, intent(in) :: n !< not used for option = 0
    integer, intent(in) :: option !< 0, 1, 2, 3, 4
    real(DP), intent(out), optional :: tottim(2) !> should be present if option=3 or 4
    integer, intent(out), optional :: nslices !> only used if option=3, still optional in that case
    
    real(DP) :: cpu,wall
    character*100 :: tmpstr

    ! no push_sub, called too frequently
    
! Check that n lies in sensible bounds

    if (n .lt. 0 .or. n .gt. mtim) then
      write(tmpstr,'(a,i6,a,i8)') 'timacc: dim mtim = ', mtim,' but input n =', n
      write(0,*) tmpstr
      stop
    end if

    if (option==0) then

! Zero out all accumulators of time and init timers
      
      acctim(:,:)=0.0d0
      tzero(:,:)=0.0d0
      ncount(:)=0

    else if (option==1) then

! Initialize timepw for n

      call timget(cpu,wall)
      tzero(1,n)=cpu
      tzero(2,n)=wall
        
    else if (option==2) then

! Accumulate time for n
      
      call timget(cpu,wall)
      acctim(1,n)=acctim(1,n)+cpu -tzero(1,n)
      acctim(2,n)=acctim(2,n)+wall-tzero(2,n)
      ncount(n)=ncount(n)+1

    else if (option==3) then

! Return accumulated time for n

      if(.not. present(tottim)) then
        write(0,*) "timacc requires tottim for option 3."
        stop
      endif

      tottim(1)=acctim(1,n)
      tottim(2)=acctim(2,n)
      if(PRESENT(nslices)) then
        nslices=ncount(n)
      end if

    else if (option==4) then

! Return elapsed time for n (do not accumulate)

      if(.not. present(tottim)) then
        write(0,*) "timacc requires tottim for option 4."
        stop
      endif

      call timget(cpu,wall)
      tottim(1)=cpu-tzero(1,n)
      tottim(2)=wall-tzero(2,n)
      
    else

      write(tmpstr,'(a,i10,a)') 'timacc: input option = ', option, 'not valid.'
      write(0,*) tmpstr
      stop

    end if

    return
  end subroutine timacc

end module timing_m
