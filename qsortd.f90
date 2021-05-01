module qsortd_mod
use kind_mod, only: dp
implicit none
private
public :: qsortd,indexx,smallest_loc,smallest,largest_loc,largest
contains
function smallest_loc(x,npick) result(indx)
! return the locations of the npick smallest elements of x(:)
real(kind=dp), intent(in) :: x(:)
integer      , intent(in) :: npick
integer                   :: indx(npick)
integer                   :: indx_all(size(x))
if (npick > size(x)) then
   write (*,*) "in smallest_loc, size(x), npick=",size(x),npick,"need npick <= size(x), STOPPING"
   error stop
end if
indx_all = indexx(x)
indx     = indx_all(:npick)
end function smallest_loc
!
function smallest(x,npick) result(y)
! return the npick smallest elements of x(:)
real(kind=dp), intent(in) :: x(:)
integer      , intent(in) :: npick
real(kind=dp)             :: y(npick)
if (npick > size(x)) then
   write (*,*) "in smallest, size(x), npick=",size(x),npick,"need npick <= size(x), STOPPING"
   error stop
end if
y = x(smallest_loc(x,npick))
end function smallest
!
function largest_loc(x,npick) result(indx)
! return the locations of the npick smallest elements of x(:)
real(kind=dp), intent(in) :: x(:)
integer      , intent(in) :: npick
integer                   :: indx(npick)
integer                   :: indx_all(size(x))
if (npick > size(x)) then
   write (*,*) "in largest_loc, size(x), npick=",size(x),npick,"need npick <= size(x), STOPPING"
   error stop
end if
indx_all = indexx(-x)
indx     = indx_all(:npick)
end function largest_loc
!
function largest(x,npick) result(y)
! return the largest elements of x(:)
real(kind=dp), intent(in) :: x(:)
integer      , intent(in) :: npick
real(kind=dp)             :: y(npick)
if (npick > size(x)) then
   write (*,*) "in largest, size(x), npick=",size(x),npick,"need npick <= size(x), STOPPING"
   error stop
end if
y = x(largest_loc(x,npick))
end function largest
!
function indexx(x) result(indx)
real(kind=dp), intent(in) :: x(:)
integer                   :: indx(size(x))
call qsortd(x,indx,size(x))
end function indexx
!
subroutine qsortd(x,ind,n)
! code converted using to_f90 by alan miller
! date: 2002-12-18  time: 11:55:47
real (dp), intent(in)  :: x(:)
integer, intent(out)   :: ind(:)
integer, intent(in)    :: n

!***************************************************************************

!                                                         robert renka
!                                                 oak ridge natl. lab.

!   this subroutine uses an order n*log(n) quick sort to sort a real (dp)
! array x into increasing order.  the algorithm is as follows.  ind is
! initialized to the ordered sequence of indices 1,...,n, and all interchanges
! are applied to ind.  x is divided into two portions by picking a central
! element t.  the first and last elements are compared with t, and
! interchanges are applied as necessary so that the three values are in
! ascending order.  interchanges are then applied so that all elements
! greater than t are in the upper portion of the array and all elements
! less than t are in the lower portion.  the upper and lower indices of one
! of the portions are saved in local arrays, and the process is repeated
! iteratively on the other portion.  when a portion is completely sorted,
! the process begins again by retrieving the indices bounding another
! unsorted portion.

! input parameters -   n - length of the array x.

!                      x - vector of length n to be sorted.

!                    ind - vector of length >= n.

! n and x are not altered by this routine.

! output parameter - ind - sequence of indices 1,...,n permuted in the same
!                          fashion as x would be.  thus, the ordering on
!                          x is defined by y(i) = x(ind(i)).

!*********************************************************************

! note -- iu and il must be dimensioned >= log(n) where log has base 2.

!*********************************************************************

integer   :: iu(21), il(21)
integer   :: m, i, j, k, l, ij, it, itt, indx
real      :: r
real (dp) :: t

! local parameters -

! iu,il =  temporary storage for the upper and lower
!            indices of portions of the array x
! m =      index for iu and il
! i,j =    lower and upper indices of a portion of x
! k,l =    indices in the range i,...,j
! ij =     randomly chosen index between i and j
! it,itt = temporary storage for interchanges in ind
! indx =   temporary index for x
! r =      pseudo random number for generating ij
! t =      central element of x

if (n <= 0) return

! initialize ind, m, i, j, and r

do  i = 1, n
  ind(i) = i
end do
m = 1
i = 1
j = n
r = .375

! top of loop

20 if (i >= j) go to 70
if (r <= .5898437) then
  r = r + .0390625
else
  r = r - .21875
end if

! initialize k

30 k = i

! select a central element of x and save it in t

ij = i + r*(j-i)
it = ind(ij)
t = x(it)

! if the first element of the array is greater than t,
!   interchange it with t

indx = ind(i)
if (x(indx) > t) then
  ind(ij) = indx
  ind(i) = it
  it = indx
  t = x(it)
end if

! initialize l

l = j

! if the last element of the array is less than t,
!   interchange it with t

indx = ind(j)
if (x(indx) >= t) go to 50
ind(ij) = indx
ind(j) = it
it = indx
t = x(it)

! if the first element of the array is greater than t,
!   interchange it with t

indx = ind(i)
if (x(indx) <= t) go to 50
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
go to 50

! interchange elements k and l

40 itt = ind(l)
ind(l) = ind(k)
ind(k) = itt

! find an element in the upper part of the array which is
!   not larger than t

50 l = l - 1
indx = ind(l)
if (x(indx) > t) go to 50

! find an element in the lower part of the array whcih is not smaller than t

60 k = k + 1
indx = ind(k)
if (x(indx) < t) go to 60

! if k <= l, interchange elements k and l

if (k <= l) go to 40

! save the upper and lower subscripts of the portion of the
!   array yet to be sorted

if (l-i > j-k) then
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  go to 80
end if

il(m) = k
iu(m) = j
j = l
m = m + 1
go to 80

! begin again on another unsorted portion of the array

70 m = m - 1
if (m == 0) return
i = il(m)
j = iu(m)

80 if (j-i >= 11) go to 30
if (i == 1) go to 20
i = i - 1

! sort elements i+1,...,j.  note that 1 <= i < j and j-i < 11.

90 i = i + 1
if (i == j) go to 70
indx = ind(i+1)
t = x(indx)
it = indx
indx = ind(i)
if (x(indx) <= t) go to 90
k = i

100 ind(k+1) = ind(k)
k = k - 1
indx = ind(k)
if (t < x(indx)) go to 100

ind(k+1) = it
go to 90
end subroutine qsortd
end module qsortd_mod