module calls

use constants
use options

contains


function flatLEblade(y)
        !returns leading and trailing edge positions for a
        !symetric blade with straight leading edge at % halfspan y
        
        !variables
        implicit none
        
        !inputs
        real :: y
        
        !output
        real, dimension(2) :: flatLEblade
        
        flatLEblade(1) = span*abs(y)*tan(sweep)
        
        flatLEblade(2) = (1.0-abs(y))*rootchord + abs(y)*tipchord
        flatLEblade(2) = flatLEblade(1) + flatLEblade(2)
end function flatLEblade


function flatTEblade(y)
        !returns leading and trailing edge positions for a
        !symmetric blade with straight trailing edge at % halfspan y
        
        !variables
        implicit none
        
        !inputs
        real :: y
        
        !outputs
        real, dimension(2) :: flatTEblade
        
        flatTEblade(2) = span*abs(y)*tan(sweep) + rootchord
        
        flatTEblade(1) =  (1.0-abs(y))*rootchord + abs(y)*tipchord
        flatTEblade(1) = -flatTEblade(1) + flatTEblade(2)
end function flatTEblade


function ellipsoid(y)
        !returns leading and trailing edge positions for a 
        !symmetric ellipsoid with separate leading / trailing
        !edge curvatures at span y
        
        !variables
        implicit none
        
        !inputs
        real :: y
        
        !outputs
        real, dimension(2) :: ellipsoid
        
        ellipsoid(2) = sqrt(1.0 - y*y)
        
        ellipsoid(1) =  leadaxis -  leadaxis*ellipsoid(2)
        ellipsoid(2) =  leadaxis + trailaxis*ellipsoid(2)
end function ellipsoid


function linspace(num,x1,x2)
        !tested
        !returns num linearly spaced points from x1 to x2 inclusive
        
        !variables
        implicit none
        
        !inputs
        real :: x1, x2
        integer :: num
        
        !outputs
        real, dimension(N1) :: linspace
        
        !internal variables
        integer :: i
        
        linspace = 0.0
        
        do 1 i = 1,num
                
                linspace(i) = x1 + (x2-x1)*(i-1.0)/(num-1.0)
        
        1 continue
end function linspace


function cosspace(num,x1,x2)
        !tested
        !returns num points with cosine spacing from x1 to x2 inclusive
        
        !variables
        implicit none
        
        !inputs
        real :: x1, x2
        integer :: num
        
        !outputs
        real, dimension(N1) :: cosspace
        
        !internal variables
        integer :: i
        
        cosspace = 0.0
        
        do 1 i = 1,num
                
                cosspace(i) = x1 + (x2-x1)*0.5*(1.0-cos(pi*(i-1.0)/(num-1.0)))
        
        1 continue
end function cosspace


function NACA2c(h,m,x)
        !tested
        !returns camber height at x% chord for NACA 2 digit camberline with peak at height h and location m

        !variables
        implicit none
        
        !inputs
        real :: h   !camber at peak (proportion of chord)
        real :: m   !point of maximum camber
        real :: x   !chord position
        
        !output
        real :: NACA2c

        !equations from NASA Technical Memorandum 4741 page 7
        !'Computer program to obtin ordinates for NACA airfoils'
        
        if (x < 0.0 .or. x > 1.0) then
                
                print *, 'NACA2c: invalid chordline coordinate'
                NACA2c = 0.0
                return
        
        else if ((m .lt. 1.0e-3) .or. (m .gt. (1.0-1.0e-3))) then
                
                print *, 'NACA2c: invalid peak position'
                NACA2c = 0.0
                return
                
        else if (x <= m) then
                
                NACA2c = h*(2.0*m*x - x*x)/(m*m)
        
        else if (x > m) then
        
                NACA2c = h*(1.0 - 2.0*m + 2.0*m*x - x*x)/((1.0 - m)*(1.0 - m))
                
        end if
end function NACA2c


subroutine NACA2c_vec(h,m,x,no,ta)
        !tested
        
        !returns normal and tangent at x% chord for NACA 2 digit camberline with peak at height h and location m
        !variables
        implicit none
        
        !inputs
        real, intent(in) :: h   !camber at peak (proportion of chord)
        real, intent(in) :: m   !point of maximum camber
        real, intent(in) :: x   !chord position

        !outputs
        real, intent(out), dimension(2) :: no   !camber normal at x
        real, intent(out), dimension(2) :: ta   !camber tangent at x
        
        !equations from NASA Technical Memorandum 4741 page 7
        !'Computer program to obtin ordinates for NACA airfoils'
        
        no(2) = 1.0
        ta(1) = 1.0
        
        if (x < 0.0 .or. x > 1.0) then
                
                print *, 'NACA2c_vec: invalid chordline coordinate'
                no = 0.0
                ta = 0.0
                return
                
        else if (x <= m) then
                
               ta(2) = 2.0*h*(m - x)/(m*m)
        
        else if (x > m) then
                
               ta(2) = 2.0*h*(m - x)/((1.0 - m)*(1.0 - m))
        
        else
                
                print *, 'invalid panel coordinates'
        
        end if
        
        no(1) = -ta(2)
        
        ta = ta/norm2(ta)
        no = no/norm2(no)
end subroutine NACA2c_vec


function vor3line(x,x1,x2)
        !tested against Luca's
        !returns the velocity induced at x by a constant unit strength vortex segment from x1 to x2
        
        !variables
        implicit none
        
        !inputs
        real, dimension(3) :: x     !point of interest
        real, dimension(3) :: x1    !start of vortex line
        real, dimension(3) :: x2    !end of vortex line
        
        !outputs
        real, dimension(3) :: u     !induced velocity at x
        real, dimension(3) :: vor3line
        
        !internal variables
        real :: rcut = 1.0e-8
        real :: vx, vy, vz, v
        real :: r1, r2, r0r1, r0r2
        
        vx =   (x(2)-x1(2))*(x(3)-x2(3)) - (x(3)-x1(3))*(x(2)-x2(2))
        
        vy = -((x(1)-x1(1))*(x(3)-x2(3)) - (x(3)-x1(3))*(x(1)-x2(1)))
        
        vz =   (x(1)-x1(1))*(x(2)-x2(2)) - (x(2)-x1(2))*(x(1)-x2(1))
        

        v = vx*vx + vy*vy + vz*vz
        
        r1 = (x(1) - x1(1))*(x(1) - x1(1)) + &
             (x(2) - x1(2))*(x(2) - x1(2)) + &
             (x(3) - x1(3))*(x(3) - x1(3))
        
        r2 = (x(1) - x2(1))*(x(1) - x2(1)) + &
             (x(2) - x2(2))*(x(2) - x2(2)) + &
             (x(3) - x2(3))*(x(3) - x2(3))
        
        r1 = sqrt(r1)
        r2 = sqrt(r2)
        
        if ((r1.lt.rcut) .or. (r2.lt.rcut) .or. (v.lt.rcut)) then
                go to 1
        end if
        
        r0r1 =  (x2(1) - x1(1))*(x(1) - x1(1)) + &
                (x2(2) - x1(2))*(x(2) - x1(2)) + &
                (x2(3) - x1(3))*(x(3) - x1(3))
        
        r0r2 =  (x2(1) - x1(1))*(x(1) - x2(1)) + &
                (x2(2) - x1(2))*(x(2) - x2(2)) + &
                (x2(3) - x1(3))*(x(3) - x2(3))
        
        r0r1 = r0r1/r1
        r0r2 = r0r2/r2
        
        v = (r0r1 - r0r2)/(4.0*pi*v)
        
        u(1) = vx*v
        u(2) = vy*v
        u(3) = vz*v
        
        vor3line = u
        
        return
        1 u = 0.0
end function vor3line


function vring(x,x1,x2,x3,x4)
        !tested against Luca's
        !returns velocity induced at x by a constant unit strength vortex ring with corners x1, x2, x3, x4
        
        !variables
        implicit none
        
        !inputs
        real, dimension(3) :: x             !point of interest
        real, dimension(3) :: x1,x2,x3,x4   !corners of ring
        
        !outputs
        real, dimension(3) :: u             !induced velocity at x
        real, dimension(3) :: vring
        
        !internal variables
        real, dimension(3) :: v
        
        v = vor3line(x,x1,x2)
        u = v
        
        v = vor3line(x,x2,x3)
        u = u + v
        
        v = vor3line(x,x3,x4)
        u = u + v
        
        v = vor3line(x,x4,x1)
        u = u + v
        
        vring = u
end function vring


function vring_trail(x,x1,x2,x3,x4)
        !tested
        !returns two velocities:
        !vring_trail(1,:) is the velocity induced at x by a constant unit strength vortex ring with corners x1, x2, x3, x4
        !vring_trail(2,:) is the velocity induced by the 'trailing' vortex segments (x2-x3 and x4-x1)

        !variables
        implicit none
        
        !inputs
        real, dimension(3) :: x             !point of interest
        real, dimension(3) :: x1,x2,x3,x4   !corners of ring
        
        !outputs
        real, dimension(3)   :: u1          !induced velocity at x
        real, dimension(3)   :: u2          !induced velocity by trailing vortices at x
        real, dimension(2,3) :: vring_trail
        
        !internal variables
        real, dimension(3) :: v
        
        u1 = 0.0
        u2 = 0.0

        v  = vor3line(x,x1,x2)
        u1 = u1 + v
        
        v  = vor3line(x,x2,x3)
        u1 = u1 + v
        u2 = u2 + v
        
        v  = vor3line(x,x3,x4)
        u1 = u1 + v
        
        v  = vor3line(x,x4,x1)
        u1 = u1 + v
        u2 = u2 + v
        
        vring_trail(1,:) = u1
        vring_trail(2,:) = u2
end function vring_trail


function indx2(i,j,base)
        !tested
        !returns consecutive integers with increasing first i, then j, where base is the maximum i value
        
        !variables
        implicit none
        
        !inputs
        integer :: i, j
        integer :: base
        
        !outputs
        integer :: indx2
        
        indx2 = (j-1)*base + i
end function indx2


function indx3(i,j,jbase,s,sbase)
        !tested
        !returns consecutive integers with increasing first i, then j, where base is the maximum i value
        
        !variables
        implicit none
        
        !inputs
        integer :: i, j, s
        integer :: jbase, sbase
        
        !outputs
        integer :: indx3
        
        indx3 = (s-1)*sbase*jbase + (j-1)*jbase + i
end function indx3


function solvecomplex(order,a,rhs)
        !tested
        !wrapper for lapack cgesv, reduces # of inputs
        
        !variables
        implicit none
        
        !function output
        complex, dimension(N2) :: solvecomplex
        
        !inputs
        complex, dimension(N2,N2) :: a     !matrix A
        complex, dimension(N2)   :: rhs   !rhs b
        integer :: order
        
        !internal variables
        integer, parameter    :: NRHS = 1               !number of right hand sides
        integer, dimension(N2) :: IPIV                   !pivot vector
        integer :: info                                 !information on solution
        
        call cgesv(order,NRHS,a,N2,IPIV,rhs,N2,info)
        
        if (info.gt.0) then
                print *, 'illegal solver value', info
        else if (info.lt.0) then
                print *, 'singular solution', info
        end if
        
        solvecomplex = rhs
end function solvecomplex


function dotcr(length,a,b)
        !tested
        !returns dot product of a complex (a) and real (b) vectors of dimension(length)
        
        !variables
        implicit none
        
        !inputs
        integer :: length
        complex, dimension(length) :: a
        real, dimension(length) :: b
        
        !outputs
        complex :: dotcr

        !internal variables
        integer :: i
        
        dotcr = 0.0
        
        do i = 1,length
                
                dotcr = dotcr + a(i)*b(i)
        
        end do
end function dotcr


function dotrc(length,a,b)
        !tested
        !returns dot product of a real (a) and complex (b) vectors of dimension(length)
        
        !variables
        implicit none
        
        !inputs
        integer :: length
        real, dimension(length) :: a
        complex, dimension(length) :: b
        
        !outputs
        complex :: dotrc

        !internal variables
        integer :: i

        dotrc = 0.0
        
        do i = 1,length
                
                dotrc = dotrc + a(i)*b(i)
        
        end do
end function dotrc


function dotcc(length,a,b)
        !tested
        !returns (non-conjugate) dot product of two complex (b) vectors of dimension(length)
        
        !variables
        implicit none
        
        !inputs
        integer :: length
        complex, dimension(length) :: a,b
        
        !outputs
        complex :: dotcc

        !internal variables
        integer :: i

        dotcc = 0.0
        
        do i = 1,length

                dotcc = dotcc + a(i)*b(i)

        end do
end function dotcc


function rvec2mat(vecN,matR,matC,vec)
        !tested
        !returns a real matrix whose rows are populated by consecutive components of vec

        !variables
        implicit none
        
        !inputs
        integer :: vecN         !number of vector components
        integer :: matR         !length of matrix rows (number of columns)
        integer :: matC         !length of matrix columns (number of rows)
        real,  dimension(vecN)   :: vec
        
        !outputs
        real, dimension(matR,matC) :: rvec2mat
        
        !internal variables
        integer :: i,j
        
        if (matR*matC .ne. vecN) then
        print *, 'rvec2mat: incompatible sizes'
        return
        end if
        
        rvec2mat = 0.0
        
        do 1 j = 1,matC
        do 1 i = 1,matR
                
                rvec2mat(i,j) = vec(indx2(i,j,matR))
        
        1 continue
end function rvec2mat


function cvec2mat(vecN,matR,matC,vec)
        !tested
        !returns a complex matrix whose rows are populated by consecutive components of vec

        !variables
        implicit none
        
        !inputs
        integer :: vecN         !number of vector components
        integer :: matR         !length of matrix rows (number of columns)
        integer :: matC         !length of matrix columns (number of rows)
        complex, dimension(vecN)   :: vec
        
        !outputs
        complex, dimension(matR,matC) :: cvec2mat
        
        !internal variables
        integer :: i

        if (matR*matC .ne. vecN) then
        print *, 'cvec2mat: incompatible sizes'
        return
        end if
        
        cvec2mat = 0.0
        
        do 1 i = 1,matR
                
                cvec2mat(i,1:matC) = vec((i-1)*matC+1:i*matC)
        
        1 continue
end function cvec2mat


function rmat2vec(vecN,matR,matC,mat)
        !tested
        !returns a real vector populated with consectutive rows of matrix mat

        !variables
        implicit none
        
        !inputs
        integer :: vecN         !number of vector components
        integer :: matR         !length of matrix rows (number of columns)
        integer :: matC         !length of matrix columns (number of rows)
        real, dimension(matR,matC) :: mat
        
        !outputs
        real, dimension(vecN)   :: rmat2vec
        
        !internal variables
        integer :: i
        
        if (matR*matC .ne. vecN) then
        print *, 'rmat2vec: incompatible sizes'
        return
        end if
        
        rmat2vec = 0.0
        
        do 1 i = 1,matR
                
                rmat2vec((i-1)*matC+1:i*matC) = mat(i,1:matC)
        
        1 continue
end function rmat2vec


function cmat2vec(vecN,matR,matC,mat)
        !tested
        !returns a complex vector populated with consectutive rows of matrix mat

        !variables
        implicit none
        
        !inputs
        integer :: vecN         !number of vector components
        integer :: matR         !length of matrix rows (number of columns)
        integer :: matC         !length of matrix columns (number of rows)
        complex, dimension(matR,matC) :: mat
        
        !outputs
        complex,  dimension(vecN)   :: cmat2vec
        
        !internal variables
        integer :: i
        
        if (matR*matC .ne. vecN) then
        print *, 'cmat2vec: incompatible sizes'
        return
        end if
        
        cmat2vec = 0.0
        
        do 1 i = 1,matR
                
                cmat2vec((i-1)*matC+1:i*matC) = mat(i,1:matC)
        
        1 continue
end function cmat2vec


function norm2c(length,a)
        !tested
        !returns the (real) norm2 length of a complex vector
        
        !variables
        implicit none
        
        !inputs
        integer :: length
        complex, dimension(length) :: a
        
        !outputs
        real :: norm2c
        
        !internal variables
        integer :: i
        
        norm2c = 0.0
        
        do i = 1,length
                
                norm2c = norm2c + real(a(i)*conjg(a(i)))
        
        end do
        
        norm2c = sqrt(norm2c)
end function norm2c


function rcross(a,b)
        !cross product of two real 3 dimensional vectors
        
        !variables
        implicit none
        
        !inputs
        real, dimension(3) :: a, b
        
        !outputs
        real, dimension(3) :: rcross
        
        rcross = 0.0
        
        rcross(1) = a(2)*b(3) - a(3)*b(2)
        rcross(2) = a(3)*b(1) - a(1)*b(3)
        rcross(3) = a(1)*b(2) - a(2)*b(1)
end function rcross
        

function ccross(a,b)
        !cross product of two complex 3 dimensional vectors
        
        !variables
        implicit none
        
        !inputs
        complex, dimension(3) :: a, b
        
        !outputs
        complex, dimension(3) :: ccross
        
        ccross = 0.0
        
        ccross(1) = a(2)*b(3) - a(3)*b(2)
        ccross(2) = a(3)*b(1) - a(1)*b(3)
        ccross(3) = a(1)*b(2) - a(2)*b(1)
end function ccross
        

end module calls
