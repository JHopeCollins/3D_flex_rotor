module constants

!variables
implicit none

!array sizes
integer, parameter :: N1 = 128
integer, parameter :: N2 = 1700

!numeric constants
real, parameter :: pi  = 4.0*atan(1.0)
real, parameter :: pi2 = 8.0*atan(1.0)
real, parameter :: d2r = pi/180.0
real, parameter :: r2d = 180.0/pi

complex, parameter :: imag = (0.0 , 1.0)

!files
integer, parameter :: file1 = 99
integer, parameter :: file2 = 98
integer, parameter :: file3 = 97
integer, parameter :: file4 = 96
integer, parameter :: file5 = 95
integer, parameter :: file6 = 94

end module constants
