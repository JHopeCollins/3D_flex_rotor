program test

use constants
use calls

!variables
implicit none


!call vortex_tests()

!call NACA2_tests()

!call indx_tests()

!call norm2c_tests()

!call space_tests()

call geo_tests()


end program

!-------------------------------------------------------

subroutine geo_tests()
        use constants
        use calls
        
        !variables
        implicit none
        
        real, dimension(N1,N1,2) :: x
        real :: x1, x2
        integer :: num
        integer :: i,j
        
        num =  15
        x1  = -1.0
        x2  =  1.0
        
        x(1,:,2) = cosspace(num,0.95*x1,0.95*x2)
        do 8 j = 1,num
                x(2:num,j,2) = x(1,j,2)
        8 continue
        
        do 1 j = 1,num
                
                x(1:num:num-1,j,1) = flatLEblade(x(1,j,2))
                x( :         ,j,1) = linspace(num, x(1,j,1), x(num,j,1))
                
                do 1 i = 1,num
                write(140,*) x(i,j,:)
        1 continue

        do 2 j = 1, num
                
                x(1:num:num-1,j,1) = flatTEblade(x(1,j,2))
                x( :         ,j,1) = linspace(num, x(1,j,1), x(num,j,1))
                
                do 2 i = 1,num
                write(141,*) x(i,j,:)
        2 continue
        
        do 3 j = 1, num
                
                x(1:num:num-1,j,1) = ellipsoid(x(1,j,2))
                x( :         ,j,1) = cosspace(num, x(1,j,1), x(num,j,1))
                
                do 3 i = 1,num
                write(142,*) x(i,j,:)
        3 continue
end subroutine geo_tests


subroutine space_tests()
        use constants
        use calls
        
        !variables
        implicit none
        integer :: num = 20
        real    :: x1 = 0.0
        real    :: x2 = 1.0
        real, dimension(N1) :: x
        integer :: i
        
        x = linspace(num,x1,x2)
        
        do 1 i = 1,num
                write(100,*) x(i), 0.0
        1 continue
        
        x = cosspace(num,x1,x2)
        
        do 2 i = 1,num
                write(101,*) x(i), 0.0
        2 continue
end subroutine space_tests
        

subroutine vortex_tests()
        use constants
        use calls
        
        !variables
        implicit none

        real, dimension(3) :: x, x1, x2, x3, x4, u
        real, dimension(3) :: xa, x1a, x2a, x3a, x4a
        real, dimension(3) :: xb, x1b, x2b, x3b, x4b
        real, dimension(2,3) :: v
        
        integer :: i,j
        
        print *, 'vortex tests'

        x   = [30.0, 30.0, 0.0]

        x1a = [-1.0, -1.0, 0.0]
        x2a = [-1.0,  1.0, 0.0]
        x3a = [ 1.0,  1.0, 0.0]
        x4a = [ 1.0, -1.0, 0.0]

        x1b = [-2.0, -2.0, 0.0]
        x2b = [-2.0,  2.0, 0.0]
        x3b = [ 2.0,  2.0, 0.0]
        x4b = [ 2.0, -2.0, 0.0]
        
        print *, 'small', vring(x,x1a,x2a,x3a,x4a)
        print *, 'large', vring(x,x1b,x2b,x3b,x4b)
        

                do 1 i = 1, 15
                do 1 j = 1, 15
                        
                        x(1) = -4.0 + 8.0*(i-1.0)/15.0
                        x(2) = -4.0 + 8.0*(j-1.0)/15.0
                        
                        u = vring(x,x1,x2,x3,x4)
                        
                        write (120,*) x(1), x(2), u(3)
                1 continue
                
                
        
                u = vor3line(x,x1,x2)
                !print *, 'vor3line', u
                
                call vortex1(x,x1,x2,u)
                !print *, 'vortex1 ', u
                
                
                u = vring(x,x1,x2,x3,x4)
                !print *, 'vring   ', u
                
                v = vring_trail(x,x1,x2,x3,x4)
                !print *, 'trail 1 ', v(1,:)
                
                call ring(x,x1,x2,x3,x4,u)
                !print *, 'ring    ', u
                
                
                !print *, 'line12  ', vor3line(x,x1,x2)
                !print *, 'line23  ', vor3line(x,x2,x3)
                !print *, 'line34  ', vor3line(x,x3,x4)
                !print *, 'line41  ', vor3line(x,x4,x1)
                !
                !        print *, ''
                !
                !v = vring_trail(x,x1,x2,x3,x4)
                !
                !print *, 'trail 2 ', v(2,:)
                !print *, 'trail t2', vring(x,x1,x2,x3,x4) - vor3line(x,x1,x2) - vor3line(x,x3,x4)
                !print *, 'trail t2', vor3line(x,x2,x3) + vor3line(x,x4,x1)
end subroutine vortex_tests


subroutine NACA2_tests()
        use constants
        use calls
        
        !variables
        implicit none
        
        real :: height, peak
        real :: x, dx, y
        real, dimension(2) :: no     , ta
        real, dimension(2) :: no_test, ta_test
        
        print *, 'NACA2 tests'
        
        height = 0.20
        peak   = 0.00
        
        x  = 0.6
        dx = 1.0e-5
        
        if ((peak.lt.1.0e-3) .or. (peak.gt.(1.0-1.0e-3))) then
                print *, 'peak'
        end if

        print *, 'height           ', height
        print *, 'peak             ', peak
        print *, ' x for no/ta test',  x
        print *, 'dx for no/ta test', dx
        
                print *, '1'
        
        print *, height, peak
        print *, 'camber @ le      '
        y = NACA2c(height,peak, 1.0)
        print *, y
                print *, '2'
        print *, 'camber @ te      ', NACA2c(height,peak, 1.0)
        print *, 'camber @ peak    ', NACA2c(height,peak,peak)
        print *, 'camber @ midchord', NACA2c(height,peak, 0.5)
        
                print *, ''

        call NACA2c_vec(height,peak,peak,no,ta)
        
        print *, 'normal  @ peak', no
        print *, 'tangent @ peak', ta
        
                print *, ''
        
        call NACA2c_vec(height,peak,x,no,ta)

        ta_test(1) =  dx
        ta_test(2) =  NACA2c(height,peak,x+dx) - NACA2c(height,peak,x)
        ta_test    =  ta_test/norm2(ta_test)
        no_test(1) = -ta_test(2)
        no_test(2) =  ta_test(1)
        
        print *, 'tangent @ x     ', ta
        print *, 'tangent test @ x', ta_test
        print *, 'normal  @ x     ', no
        print *, 'normal  test @ x', no_test
end subroutine NACA2_tests


subroutine indx_tests()
        use constants
        use calls
        
        !variables
        implicit none
        
        integer :: i,j,base
        
        print *, 'indx tests'
        
        base = 5
        print *, 'base', base
        
                print *, ''
        
        do 1 j = 1,3
        do 1 i = 1,base
                
                print *, 'i,j, indx', i, j, indx(i,j,base)
                
        1 continue
end subroutine indx_tests


subroutine norm2c_tests()
        use constants
        use calls
        
        !variables
        implicit none
        
        complex, dimension(2) :: a
        
        a(1) = (1.0,2.0)
        a(2) = (3.0,4.0)
        
        print *, 'norm2c', norm2c(2,a)
        print *, 'test  ', sqrt(real(a(1))**2.0 + aimag(a(1))**2.0 + real(a(2))**2.0 + aimag(a(2))**2.0)
end subroutine norm2c_tests


!-------------------------------------------------------

!Luca's functions to verify against

subroutine vortex1(x,x1,x2,u)
        implicit none
        
        real              pi
        parameter         (pi=3.1415926585)
  
        real              x,x1,x2,u
        dimension         x(3), x1(3), x2(3), u(3)
  
        real              rcut,vx,vy,vz,v2,r1,r2,r0r1,r0r2, &
                         coef
        
        
        rcut=1.0e-10
        
        !vx=  (y   -y1)   *(z   -z2)   -(z   -z1)   *(y   -y2)
        vx=  (x(2)-x1(2))*(x(3)-x2(3))-(x(3)-x1(3))*(x(2)-x2(2))
  
        !vy=-((x   -x1)   *(z   -z2)   -(z   -z1)   *(x   -x2))
        vy=-((x(1)-x1(1))*(x(3)-x2(3))-(x(3)-x1(3))*(x(1)-x2(1)))
  
        !vz=  (x   -x1)   *(y   -y2)   -(y   -y1)   *(x   -x2)
        vz=  (x(1)-x1(1))*(x(2)-x2(2))-(x(2)-x1(2))*(x(1)-x2(1))
        
        v2=vx*vx+vy*vy+vz*vz
        !r1=sqrt((x   -x1)   *(x   -x1)   +(y   -y1)   *(y   -y1)   +(z   -z1)   *(z   -z1))
        r1=sqrt((x(1)-x1(1))*(x(1)-x1(1))+(x(2)-x1(2))*(x(2)-x1(2))+(x(3)-x1(3))*(x(3)-x1(3)))
  
        !r2=sqrt((x   -x2)   *(x   -x2)   +(y   -y2)   *(y   -y2)   +(z   -z2)   *(z   -z2))
        r2=sqrt((x(1)-x2(1))*(x(1)-x2(1))+(x(2)-x2(2))*(x(2)-x2(2))+(x(3)-x2(3))*(x(3)-x2(3)))
        
        if((r1 .lt. rcut) .or. (r2 .lt. rcut) .or. (v2 .lt. rcut) ) goto 1 
        
        !r0r1=(x2   -x1)   *(x   -x1)   +(y2   -y1)   *(y   -y1)   +(z2   -z1)   *(z   -z1)
        r0r1=(x2(1)-x1(1))*(x(1)-x1(1))+(x2(2)-x1(2))*(x(2)-x1(2))+(x2(3)-x1(3))*(x(3)-x1(3))
  
        !r0r2=(x2   -x1   )*(x   -x2)   +(y2   -y1)   *(y   -y2)   +(z2   -z1)   *(z   -z2)
        r0r2=(x2(1)-x1(1))*(x(1)-x2(1))+(x2(2)-x1(2))*(x(2)-x2(2))+(x2(3)-x1(3))*(x(3)-x2(3))
  
        coef=(r0r1/r1-r0r2/r2)/(4.0*pi*v2)
  
        u(1)=vx*coef
        u(2)=vy*coef
        u(3)=vz*coef
        
        goto 2
        
        
        1 u=0.
        
        2 continue
  
        return
end 


subroutine ring(x,x1,x2,x3,x4,uf)
        implicit none
       
        real      x,x1,x2,x3,x4,uf,u
        dimension x(3),x1(3),x2(3),x3(3),x4(3),uf(3),u(3)
       
       
        call vortex1( x, x1, x2, u)
       
        uf(1)=u(1)
        uf(2)=u(2)
        uf(3)=u(3)
       
        call vortex1( x, x2, x3, u)
       
        uf(1)=uf(1)+u(1)
        uf(2)=uf(2)+u(2)
        uf(3)=uf(3)+u(3)
       
        call vortex1( x, x3, x4, u)
       
        uf(1)=uf(1)+u(1)
        uf(2)=uf(2)+u(2)
        uf(3)=uf(3)+u(3)
       
        call vortex1( x, x4, x1, u)
       
        uf(1)=uf(1)+u(1)
        uf(2)=uf(2)+u(2)
        uf(3)=uf(3)+u(3)
end 
