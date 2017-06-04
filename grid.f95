!construct grid point (panel corners) array, and panel geometry
subroutine grid(x,xc,xv,xn,xtc,xts)

use constants
use calls
use options

!variables
implicit none

!outputs
real, intent(out), dimension(Nb,N1,N1,3) :: x        !panel corners
real, intent(out), dimension(Nb,N1,N1,3) :: xc       !collocation points
real, intent(out), dimension(Nb,N1,N1,3) :: xv       !vortex ring corners
real, intent(out), dimension(Nb,N1,N1,3) :: xn       !panel normal vectors
real, intent(out), dimension(Nb,N1,N1,3) :: xtc      !panel tangent vectors chordwise
real, intent(out), dimension(Nb,N1,N1,3) :: xts      !panel tangent vectors spanwise

!internal variables
integer :: i,j, s
real :: xspan   !local spanwise coordinate
real :: xchord  !local chordwise coordinate
real :: chord   !local chord
real :: theta   !angle of wake point
real :: dtheta  !rotor rotation in one timestep
real :: radius  !radius of wake point
real :: angle   !angle to rotate blade geometry

real, dimension(3,3) :: R

real, dimension(3) :: Ak, Bk    !panel diagonals

!initialise
x  = 0.0
xc = 0.0
xv = 0.0
xn = 0.0
Ak = 0.0
Bk = 0.0
xchord = 0.0
xspan  = 0.0
chord  = 0.0
dtheta = dt*Romega

!generate single blade geometry

!panel corners
 x(1,1,:,2) = linspace(Ns+1, hub, hub + clip*AR)
!x(1,1,:,2) = cosspace(Ns+1, hub, hub + clip*AR)

do 1 j = 1,Ns+1
        
        x(1,2:Nc+1,j,2)    = x(1,1,j,2)
        
        x(1,1:Nc+1:Nc,j,1) = flatLEblade((x(1,1,j,2)-hub)/AR)
        !x(1,1:Nc+1:Nc,j,1) = flatTEblade((x(1,1,j,2)-hub)/AR)
        !x(1,1:Nc+1:Nc,j,1) = ellipsoid(  (x(1,1,j,2)-hub)/AR)
        
         x(1,:,j,1) = linspace(Nc+1, x(1,1,j,1), x(1,Nc+1,j,1))
        !x(1,:,j,1) = cosspace(Nc+1, x(1,1,j,1), x(1,Nc+1,j,1))
        
        chord = x(1,Nc+1,j,1) - x(1,1,j,1)
        
        do 1 i = 1,Nc+1
                
                xchord = (x(1,i   ,j,1) - x(1,1,j,1)) / (1.001*chord)
                
                 x(1,i,j,3) = NACA2c(h,m,xchord)
                 x(1,i,j,3) = chord*x(1,i,j,3)
                !x(1,i,j,3) = 0.0

1 continue

!collocation points at panel centre span, 3/4 chord
do 3 j = 1,Ns
do 3 i = 1,Nc
        
        xc(1,i,j,1) = 0.125*( x(1,i  ,j,1) + x(1,i  ,j+1,1) ) + &
                      0.375*( x(1,i+1,j,1) + x(1,i+1,j+1,1) )
        
        xc(1,i,j,2) = 0.25 *( x(1,i  ,j,2) + x(1,i+1,j  ,2)   + &
                              x(1,i,j+1,2) + x(1,i+1,j+1,2) )
        
        xc(1,i,j,3) = 0.125*( x(1,i  ,j,3) + x(1,i  ,j+1,3) ) + &
                      0.375*( x(1,i+1,j,3) + x(1,i+1,j+1,3) )

3 continue

blade = xc(1,:,:,:)

!bound vortex ring corners. leading vortex line along panel 1/4 chord
do 4 j = 1,Ns+1
do 4 i = 1,Nc
        
        xv(1,i,j,1) = 0.75*x(1,i,j,1) + 0.25*x(1,i+1,j,1)
        xv(1,i,j,2) = 0.75*x(1,i,j,2) + 0.25*x(1,i+1,j,2)
        xv(1,i,j,3) = 0.75*x(1,i,j,3) + 0.25*x(1,i+1,j,3)
        
        
        !xv(1,i,j,1:3) = 0.75*x(1,i,j,1:3) + 0.25*x(1,i+1,j,1:3)
        
4 continue


!pitch rotation
R(1,1:3) = [ cos(pitch) , 0.0 , sin(pitch)]
R(2,1:3) = [    0.0     , 1.0 ,    0.0    ]
R(3,1:3) = [-sin(pitch) , 0.0 , cos(pitch)]

do 10 j = 1,Ns+1
        
        do 10 i = 1,Nc+1
        
        x( 1,i,j,1:3) = matmul(R, x(1,i,j,1:3))
        xc(1,i,j,1:3) = matmul(R,xc(1,i,j,1:3))
        xv(1,i,j,1:3) = matmul(R,xv(1,i,j,1:3))

10 continue


!generate all blades with identical geometry
!rotation to space blades equally around axis
do 20 s = Nb,1,-1

        angle = pi/2.0 + (s-1.0)*2.0*pi/Nb
        
        R(1,1:3) = [ 1.0 ,    0.0    ,     0.0    ]
        R(2,1:3) = [ 0.0 , cos(angle), -sin(angle)]
        R(3,1:3) = [ 0.0 , sin(angle),  cos(angle)]

        do 20 j = 1,Ns+1
        do 20 i = 1,Nc+1
                
                x( s,i,j,1:3) = matmul(R, x(1,i,j,1:3))
                xc(s,i,j,1:3) = matmul(R,xc(1,i,j,1:3))
                xv(s,i,j,1:3) = matmul(R,xv(1,i,j,1:3))
                
20 continue

!panel normal and tangent vectors
do 2 s = 1,Nb
do 2 j = 1,Ns
do 2 i = 1,Nc
        
        !panel diagonals
         Ak(:) = x(s,i+1,j+1,:) - x(s,i  ,j,:)
         Bk(:) = x(s,i  ,j+1,:) - x(s,i+1,j,:)
        
        !normals by cross product
        xn(s,i,j,:) = rcross(Ak,Bk)
        !xn(s,i,j,1) = Ak(2)*Bk(3) - Ak(3)*Bk(2)
        !xn(s,i,j,2) = Ak(3)*Bk(1) - Ak(1)*Bk(3)
        !xn(s,i,j,3) = Ak(1)*Bk(2) - Ak(2)*Bk(1)
        
        !chordwise tangent
        xtc(s,i,j,:) =  0.5*(x(s,i+1,j  ,:) - x(s,i,j  ,:)  + &
                             x(s,i+1,j+1,:) - x(s,i,j+1,:))
        
        !spanwise tangent
        xts(s,i,j,:) = 0.25*(x(s,i  ,j+1,:) - x(s,i  ,j,:)) + &
                       0.75*(x(s,i+1,j+1,:) - x(s,i+1,j,:))

        xn( s,i,j,:) = xn( s,i,j,:)/norm2(xn( s,i,j,:))
        xtc(s,i,j,:) = xtc(s,i,j,:)/norm2(xtc(s,i,j,:))
        xts(s,i,j,:) = xts(s,i,j,:)/norm2(xts(s,i,j,:))

2 continue


!shed vortex ring corners. vortex lines convected with mean flow
do 5 s = 1,Nb
do 5 j = 1,Ns+1
        
        theta = atan(xv(s,Nc,j,3)/xv(s,Nc,j,2))
        if (xv(s,Nc,j,2).lt.0.0) then
                theta = theta + pi
        end if
        
        radius = xv(s,Nc,j,2)*xv(s,Nc,j,2) + xv(s,Nc,j,3)*xv(s,Nc,j,3)
        radius = sqrt(radius)
        
        do 5 i = Nc+1,Nw
                
                theta = theta + dtheta
                
                xv(s,i,j,:) = xv(s,i-1,j,:) + dt*uinf(:)
                
                xv(s,i,j,2) = radius*cos(theta)
                xv(s,i,j,3) = radius*sin(theta)
                
5 continue


do 99 j = 1,Ns+1
        do 98 i = 1,Nc+1
        write(120,*)  x(1,i,j,:), &
                      x(2,i,j,:), &
                      x(3,i,j,:)

        write(121,*) xc(1,i,j,:), &
                     xc(2,i,j,:), &
                     xc(3,i,j,:)

        if (i.eq.j) go to 98
                
        write(122,*) xv(1,i,j,:), &
                     xv(2,i,j,:), &
                     xv(3,i,j,:)
                
        98 continue
        write (120,*) ''
        write (121,*) ''
        write (122,*) ''
99 continue

do 95 j = 1,Ns+1
        do 94 i = Nc,Nw
                write(123,*) xv(1,i,j,:), xv(2,i,j,:), xv(3,i,j,:) 
        94 continue
        write(123,*) ''
95 continue

end subroutine grid


!---------------------------------------------------------------------------


!change grid points (panel corners) due to blade flexing
subroutine flex_grid(x, xn, dx, dxc, dxv, dxn)

use constants
use calls
use options

!variables
implicit none

!inputs
real,    intent(in),  dimension(Nb,N1,N1,3) :: x         !panel corners
real,    intent(in),  dimension(Nb,N1,N1,3) :: xn        !panel corners

!outputs
complex, intent(out), dimension(Nb,N1,N1,3) :: dx        !panel corners
complex, intent(out), dimension(Nb,N1,N1,3) :: dxc       !collocation points
complex, intent(out), dimension(Nb,N1,N1,3) :: dxv       !vortex ring corners
complex, intent(out), dimension(Nb,N1,N1,3) :: dxn       !panel normal vectors

!internal variables
integer :: i, j, s
real :: phase
real, dimension(3) ::  Ak,  Bk    !panel diagonals
complex, dimension(3) :: dAk, dBk    !panel diagonals
real, dimension(2) :: no,ta

!initialise
        no  = 0.0
        ta  = 0.0
        dx  = 0.0
        dxc = 0.0
        dxv = 0.0
        dxn = 0.0
        Ak  = 0.0
        Bk  = 0.0
        dAk = 0.0
        dBk = 0.0


do 1 s = 1,Nb
        
        !corners
        
                !phase = 0.0
                phase = x(s,1,1,1)
                phase = omega*phase
                dx(s,1,1,:) = xn(s,1,1,:)*exp(imag*phase)
                
                !phase = 0.0
                phase = x(s,1,Ns+1,1)
                phase = omega*phase
                dx(s,1,Ns+1,:) = xn(s,1,Ns,:)*exp(imag*phase)
        
                !phase = norm2(x(s,Nc+1,1,:)-x(s,1,1,:))
                phase =        x(s,Nc+1,1,1)
                phase = omega*phase
                dx(s,Nc+1,1,:) = xn(s,Nc,1,:)*exp(imag*phase)
        
                !phase = norm2(x(s,Nc+1,Ns+1,:)-x(s,1,Ns+1,:))
                phase =        x(s,Nc+1,Ns+1,1)
                phase = omega*phase
                dx(s,Nc+1,Ns+1,:) = xn(s,Nc,Ns,:)*exp(imag*phase)
        
        do 11 j = 2,Ns
                
                !leading edge
                !phase = 0.0
                phase = x(s,1,j,1)
                phase = omega*phase
                dx(s,1,j,:) = 0.5*( xn(s,1,j,:) + xn(s,1,j-1,:))
                dx(s,1,j,:) =       dx(s,1,j,:)*exp(imag*phase)
                
                !trailing edge
                !phase = 0.0
                phase = x(s,Nc+1,j,1)
                phase = omega*phase
                dx(s,Nc+1,j,:) = 0.5*( xn(s,Nc+1,j,:) + xn(s,Nc+1,j-1,:))
                dx(s,Nc+1,j,:) =       dx(s,Nc+1,j,:)*exp(imag*phase)
                
        11 continue
        
        do 1 i = 2,Nc
                
                !root edge
                !phase = norm2(x(s,i,1,:)-x(s,1,1,:))
                phase =        x(s,i,1,1)
                phase = omega*phase
                
                dx(s,i,1,:) = 0.5*( xn(s,i,1,:) + xn(s,i-1,1,:))
                dx(s,i,1,:) =       dx(s,i,1,:)*exp(imag*phase)
        
                !tip edge
                !phase = norm2(x(s,i,Ns+1,:)-x(s,1,Ns+1,:))
                phase =        x(s,i,Ns+1,1)
                phase = omega*phase
                
                dx(s,i,Ns+1,:) = 0.5*( xn(s,i,Ns+1,:) + xn(s,i-1,Ns+1,:))
                dx(s,i,Ns+1,:) =       dx(s,i,Ns+1,:)*exp(imag*phase)
                
                do 1 j = 2,Ns
                        
                        !phase = norm2(x(s,i,j,:) - x(s,1,j,:))
                        phase = x(s,i,j,1)
                        phase = omega*phase
                        
                        dx(s,i,j,:) = 0.25*(xn(s,i-1,j-1,:) + &
                                            xn(s,i  ,j-1,:) + &
                                            xn(s,i-1,j  ,:) + &
                                            xn(s,i  ,j  ,:))
                        
                        dx(s,i,j,:) = dx(s,i,j,:)*exp(imag*phase)
        
1 continue


!change in panel normals
! d(a x b) = (da x b) + (a x db)
do 2 s = 1,Nb
do 2 j = 1,Ns
do 2 i = 1,Nc
        
        !panel diagonals
        Ak(1) = x(s,i+1,j+1,1) - x(s,i,j,1)
        Ak(2) = x(s,i+1,j+1,2) - x(s,i,j,2)
        Ak(3) = x(s,i+1,j+1,3) - x(s,i,j,3)
        
        Bk(1) = x(s,i,j+1,1) - x(s,i+1,j,1)
        Bk(2) = x(s,i,j+1,2) - x(s,i+1,j,2)
        Bk(3) = x(s,i,j+1,3) - x(s,i+1,j,3)
        
        !panel diagonals
        dAk(1) = dx(s,i+1,j+1,1) - dx(s,i,j,1)
        dAk(2) = dx(s,i+1,j+1,2) - dx(s,i,j,2)
        dAk(3) = dx(s,i+1,j+1,3) - dx(s,i,j,3)
        
        dBk(1) = dx(s,i,j+1,1) - dx(s,i+1,j,1)
        dBk(2) = dx(s,i,j+1,2) - dx(s,i+1,j,2)
        dBk(3) = dx(s,i,j+1,3) - dx(s,i+1,j,3)
        
        !cross product diags
        dxn(s,i,j,1) = dAk(2)*Bk(3)-dAk(3)*Bk(2) + Ak(2)*dBk(3)-Ak(3)*dBk(2)
        dxn(s,i,j,2) = dAk(3)*Bk(1)-dAk(1)*Bk(3) + Ak(3)*dBk(1)-Ak(1)*dBk(3)
        dxn(s,i,j,3) = dAk(1)*Bk(2)-dAk(2)*Bk(1) + Ak(1)*dBk(2)-Ak(2)*dBk(1)
        
        dxn(s,i,j,:) = dxn(s,i,j,:)/norm2c(3,dxn(s,i,j,:))

2 continue


!change in collocation points at panel centre span, 3/4 chord
do 3 s = 1,Nb
do 3 j = 1,Ns
do 3 i = 1,Nc
        
        dxc(s,i,j,:) = 0.125*dx(s,i,j  ,:) + 0.375*dx(s,i+1,j  ,:) + &
                       0.125*dx(s,i,j+1,:) + 0.375*dx(s,i+1,j+1,:)
        
3 continue


!change in bound vortex ring corners. leading vortex line along panel 1/4 chord
do 4 s = 1,Nb
do 4 j = 1,Ns+1
do 4 i = 1,Nc
        
        dxv(s,i,j,:) = 0.75*dx(s,i,j,:) + 0.25*dx(s,i+1,j,:)

4 continue


!change in shed vortex ring corners. vortex lines convected with mean flow
do 5 s = 1,Nb
do 5 j = 1,Ns+1
do 5 i = Nc+1,Nw
        
        dxv(s,i,j,1) = 0.0
        dxv(s,i,j,2) = 0.0
        dxv(s,i,j,3) = 0.0

5 continue


end subroutine flex_grid
