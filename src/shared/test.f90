program testsac

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points
!
! returns: sorted indexing array (ibool),  reordering array (iglob) &
! number of global points (nglob)

  implicit none

  integer :: npointot

  double precision, dimension(:), allocatable :: x1, y1, z1
  integer, dimension(:), allocatable :: ibool1, iglob1, locval1, ninseg1
  logical, dimension(:), allocatable :: ifseg1
  integer :: nglob1

  double precision, dimension(:), allocatable :: x2, y2, z2
  integer, dimension(:), allocatable :: ibool2, iglob2, locval2
  integer :: nglob2

  double precision, parameter :: xtol = 1e-13

  character(len=256) :: arg
  integer :: i
  double precision :: start_time, end_time

  if (command_argument_count() < 1) then
    stop 'Missing point count.'
  endif

  call get_command_argument(1, arg)
  read(arg,*) npointot

  allocate(x1(npointot),y1(npointot),z1(npointot))
  allocate(ibool1(npointot),iglob1(npointot),locval1(npointot),ninseg1(npointot))
  allocate(ifseg1(npointot))

  allocate(x2(npointot),y2(npointot),z2(npointot))
  allocate(ibool2(npointot),iglob2(npointot),locval2(npointot))

  call random_number(x1(1:npointot/2))
  call random_number(y1(1:npointot/4))
  call random_number(z1(1:npointot/2))
  x1(npointot/2+1:npointot) = x1(1:npointot/2)
  y1(npointot/4+1:npointot/4*2) = y1(1:npointot/4)
  y1(npointot/4*2+1:npointot/4*3) = y1(1:npointot/4:-1)
  y1(npointot/4*3+1:npointot) = y1(1:npointot/4)
  z1(npointot/2+1:npointot) = z1(1:npointot/2)
  ibool1(:) = 0
  iglob1(:) = 0
  locval1(:) = 0
  ninseg1(:) = 0

  x2(:) = x1(:)
  y2(:) = y1(:)
  z2(:) = z1(:)
  ibool2(:) = 0
  iglob2(:) = 0
  locval2(:) = 0

  call cpu_time(start_time)
  call sort_array_coordinates(npointot,x1,y1,z1,ibool1,iglob1,locval1,ifseg1, &
                              nglob1,ninseg1,xtol)
  call cpu_time(end_time)
  print *, 'Test 1', end_time - start_time

  call cpu_time(start_time)
  call sort_array_coordinates2(npointot,x2,y2,z2,ibool2,iglob2,locval2, &
                               nglob2,xtol)
  call cpu_time(end_time)
  print *, 'Test 2', end_time - start_time

  print *, 'x', all(x1 == x2)
  print *, 'y', all(y1 == y2)
  print *, 'z', all(z1 == z2)
  print *, 'ibool', all(ibool1 == ibool2)
  print *, 'nglob', nglob1 == nglob2
  print *, 'iglob', all(iglob1 == iglob2)
  print *, 'locval', all(locval1 == locval2)

  open(unit=1,file='sac1')
  write(1,*) nglob1
  open(unit=2,file='sac2')
  write(2,*) nglob2
  if (command_argument_count() > 1) then
    do i=1,npointot
      write(1,'(3(g30.24,1x))') &
        x1(i),y1(i),z1(i)
      write(2,'(3(g30.24,1x))') &
        x2(i),y2(i),z2(i)
    enddo
  else
    do i=1,npointot
      write(1,'(3(g16.10,1x),3(i10,1x),l1)') &
        x1(i),y1(i),z1(i),ibool1(i),iglob1(i),locval1(i)!,ifseg1(i)
      write(2,'(3(g16.10,1x),3(i10,1x),l1)') &
        x2(i),y2(i),z2(i),ibool2(i),iglob2(i),locval2(i)
    enddo
  endif
  close(1)
  close(2)

  deallocate(x1,y1,z1)
  deallocate(ibool1,iglob1,locval1,ninseg1)
  deallocate(ifseg1)

  deallocate(x2,y2,z2)
  deallocate(ibool2,iglob2,locval2)

contains

  subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t
  
    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  end subroutine init_random_seed

  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    use iso_fortran_env, only: int64
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg

end program testsac
