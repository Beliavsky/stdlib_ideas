program xqsortd
use qsortd_mod, only: qsortd,smallest_loc,smallest,largest_loc,largest,indexx
implicit none
integer, parameter :: dp = kind(1.0d0), n = 20, npick = 5
real(kind=dp)      :: x(n)
character (len=*), parameter :: fmt_cr = "(a15,1000f8.5)"
call random_number(x)
print*,"n =",n
write (*,fmt_cr) "raw",x
write (*,fmt_cr) "ascending",x(indexx(x))
write (*,fmt_cr) "descending",x(indexx(-x))
write (*,fmt_cr) "smallest",x(smallest_loc(x,npick))
write (*,fmt_cr) "smallest",smallest(x,npick)
write (*,fmt_cr) "largest",x(largest_loc(x,npick))
write (*,fmt_cr) "largest",largest(x,npick)
end program xqsortd