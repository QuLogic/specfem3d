!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! subroutines to sort MPI buffers to assemble between chunks

  subroutine sort_array_coordinates2(npointot,x,y,z,ibool,iglob,locval, &
                                     nglob,xtol)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points
!
! returns: sorted indexing array (ibool),  reordering array (iglob) & number of global points (nglob)

  implicit none

  integer, parameter :: NDIM = 3

  integer, intent(in) :: npointot
  double precision, dimension(npointot), intent(inout) :: x, y, z
  integer, dimension(npointot), intent(inout) :: ibool
  integer, dimension(npointot), intent(out) :: iglob, locval
  integer, intent(out) :: nglob
  double precision, intent(in) :: xtol

  ! local parameters
  integer :: i
  integer :: ig

  ! establish initial pointers
  do i=1,npointot
    locval(i) = i
  enddo

  call heap_sort_multi2(npointot, x, y, z, ibool, locval, xtol)

  ! assign global node numbers (now sorted lexicographically)
  ig = 1
  iglob(locval(1)) = ig
  do i=2,npointot
    ! check for jumps in current coordinate
    if (dabs(x(i) - x(i-1)) > xtol) then
      ig = ig + 1
    elseif (dabs(y(i) - y(i-1)) > xtol) then
      ig = ig + 1
    elseif (dabs(z(i) - z(i-1)) > xtol) then
      ig = ig + 1
    endif
    iglob(locval(i)) = ig
  enddo

  nglob = ig

  end subroutine sort_array_coordinates2

! -------------------- library for sorting routine ------------------

  subroutine heap_sort_multi2(N, dx, dy, dz, ia, ib, tol)

  implicit none
  integer, intent(in) :: N
  double precision, dimension(N), intent(inout) :: dx
  double precision, dimension(N), intent(inout) :: dy
  double precision, dimension(N), intent(inout) :: dz
  integer, dimension(N), intent(inout) :: ia
  integer, dimension(N), intent(inout) :: ib
  double precision, intent(in) :: tol

  integer :: i

  ! checks if anything to do
  if (N < 2) return

  ! builds heap
  do i = N/2, 1, -1
    call heap_sort_siftdown2(i, n)
  enddo

  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    call dswap2(dx, 1, i)
    call dswap2(dy, 1, i)
    call dswap2(dz, 1, i)
    call iswap2(ia, 1, i)
    call iswap2(ib, 1, i)
    call heap_sort_siftdown2(1, i - 1)
  enddo

  contains

    subroutine dswap2(A, i, j)

    double precision, dimension(:), intent(inout) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j

    double precision :: tmp

    tmp = A(i)
    A(i) = A(j)
    A(j) = tmp

    end subroutine

    subroutine iswap2(A, i, j)

    integer, dimension(:), intent(inout) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j

    integer :: tmp

    tmp = A(i)
    A(i) = A(j)
    A(j) = tmp

    end subroutine

    subroutine heap_sort_siftdown2(start, bottom)

    integer, intent(in) :: start
    integer, intent(in) :: bottom

    integer :: i, j
    double precision :: xtmp, ytmp, ztmp
    integer :: atmp, btmp
    double precision :: xdiff, ydiff, zdiff

    i = start
    xtmp = dx(i)
    ytmp = dy(i)
    ztmp = dz(i)
    atmp = ia(i)
    btmp = ib(i)

    j = 2 * i
    do while (j <= bottom)
      ! chooses larger value first in this section
      if (j < bottom) then
        xdiff = dx(j+1) - dx(j)
        if (xdiff > tol) then
          j = j + 1
        else
          ydiff = dy(j+1) - dy(j)
          if (dabs(xdiff) <= tol .and. ydiff > tol) then
            j = j + 1
          else
            zdiff = dz(j) - dz(j+1)
            if (dabs(xdiff) <= tol .and. dabs(ydiff) <= tol .and. zdiff <= tol) then
              j = j + 1
            endif
          endif
        endif
      endif

      ! checks if section already smaller than initial value
      xdiff = xtmp - dx(j)
      if (xdiff > tol) exit
      ydiff = ytmp - dy(j)
      if (dabs(xdiff) <= tol .and. ydiff > tol) exit
      zdiff = dz(j) - ztmp
      if (dabs(xdiff) <= tol .and. dabs(ydiff) <= tol .and. zdiff > tol) exit

      dx(i) = dx(j)
      dy(i) = dy(j)
      dz(i) = dz(j)
      ia(i) = ia(j)
      ib(i) = ib(j)
      i = j
      j = 2 * i
    enddo

    dx(i) = xtmp
    dy(i) = ytmp
    dz(i) = ztmp
    ia(i) = atmp
    ib(i) = btmp

    end subroutine heap_sort_siftdown2

  end subroutine heap_sort_multi2
