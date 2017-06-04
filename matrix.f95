!subroutine to construct influence matrix and rhs of impermeability
!condition

subroutine matrix(x, xc, xv, xn, xtc, xts, a,rhs, aGc, aGs, alc, als)

use constants
use calls
use options

!variables
implicit none

!inputs
real, intent(in), dimension(Nb,N1,N1,3) :: x         !grid points
real, intent(in), dimension(Nb,N1,N1,3) :: xc        !collocation points
real, intent(in), dimension(Nb,N1,N1,3) :: xv        !vortex ring corners
real, intent(in), dimension(Nb,N1,N1,3) :: xn        !panel normals
real, intent(in), dimension(Nb,N1,N1,3) :: xtc       !panel chordwise tangent
real, intent(in), dimension(Nb,N1,N1,3) :: xts       !panel spanwise tangent

!outputs
real, intent(out), dimension(N2,N2) :: a          !influence matrix
real, intent(out), dimension(N2)   :: rhs        !rhs of linear system

real, intent(out), dimension(N2,N2) :: alc, als   !influence matrix for self induced tangent velocities
real, intent(out), dimension(N2,N2) :: aGc, aGs   !influence matrix for tangent velocities due to every other vortex ring

!internal variables
integer :: i,j, p,q, r,s
integer :: countC,countV
real, dimension(3) :: u     !velocities
real, dimension(3) :: urot  !rotor rotational velocity
real :: dx, dy              !panel length, width

!initialise
u   = 0.0
a   = 0.0
rhs = 0.0
aGc = 0.0
aGs = 0.0

!influence matrix arranged in blocks of the blade section at each spanwise location
!scan over collocation points
do 1 s = 1,Nb
do 1 j = 1,Ns
do 1 i = 1,Nc
        
        countC = indx3(i,j,Nc,s,Ns)
        
        urot = Romega*[0.0 , -xc(s,i,j,3), xc(s,i,j,2)]
        
        rhs(countC) = -dot_product(uinf , xn(s,i,j,:)) &
                      -dot_product(urot , xn(s,i,j,:))
        
        do 1 r = 1,Nb
        do 1 q = 1,Ns

                !scan over bound vortex rings
                do 2 p = 1,Nc
                        
                        countV = indx3(p,q,Nc,r,Ns)
                        
                        u = vring(xc(s,i,j,:), xv(r,p  ,q  ,:), xv(r,p  ,q+1,:), &
                                               xv(r,p+1,q+1,:), xv(r,p+1,q  ,:))
                        
                        a(countC,countV) = a(countC,countV) + dot_product(u , xn(s,i,j,:))
                        
                        aGc(countC,countV) = aGc(countC,countV) + dot_product(u, xtc(s,i,j,:))
                        aGs(countC,countV) = aGs(countC,countV) + dot_product(u, xts(s,i,j,:))
                        
                2 continue
                
                countV = indx3(Nc,q,Nc,r,Ns)

                !scan over shed vortex rings
                do 1 p = Nc+1,Nw-1
                        
                        u = vring(xc(s,i,j,:), xv(r,p  ,q  ,:), xv(r,p  ,q+1,:), &
                                               xv(r,p+1,q+1,:), xv(r,p+1,q  ,:))
                        
                        a(countC,countV) = a(countC,countV) + dot_product(u , xn(s,i,j,:))
        
                        aGc(countC,countV) = aGc(countC,countV) + dot_product(u, xtc(s,i,j,:))
                        aGs(countC,countV) = aGs(countC,countV) + dot_product(u, xts(s,i,j,:))
                        
1 continue

do 3 s = 1,Nb
do 3 j = 1,Ns
do 3 i = 1,Nc
        
        countC = indx3(i,j,Nc,s,Ns)
        
        dy = norm2(0.5*(x(s,i  ,j+1,:) - x(s,i  ,j,:) + &
                        x(s,i+1,j+1,:) - x(s,i+1,j,:)))
        
        dx = norm2(0.5*(x(s,i+1,j  ,:) - x(s,i,j  ,:) + &
                        x(s,i+1,j+1,:) - x(s,i,j+1,:)))
        
        u = vring(xc(s,i,j,:), xv(s,i  ,j  ,:), xv(s,i  ,j+1,:), &
                               xv(s,i+1,j+1,:), xv(s,i+1,j  ,:))
        
        aGc(countC,countC) = aGc(countC,countC) - dot_product(u, xtc(s,i,j,:))
        aGs(countC,countC) = aGs(countC,countC) - dot_product(u, xts(s,i,j,:))
        
        alc(countC,countC)   =  1.0/dx
        als(countC,countC)   =  1.0/dy

        if (i.ne.1) then

                alc(countC,countC-1) = -1.0/dx
        end if
        
        if (j.ne.1) then
                
                als(countC,countC-Nc) = -1.0/dy
        end if

3 continue


end subroutine matrix

!--------------------------------------------------------------------------


!influence matrix and rhs of flexible impermeability
!condition

subroutine flex_matrix(xc, xv, xn, xtc, xts, dxc, dxn, da, drhs, daGc, daGs)

use constants
use calls
use options

!variables
implicit none

!inputs
real,    intent(in), dimension(Nb,N1,N1,3) :: xc        !collocation points
real,    intent(in), dimension(Nb,N1,N1,3) :: xv        !vortex ring corners
real,    intent(in), dimension(Nb,N1,N1,3) :: xn        !panel normals
real,    intent(in), dimension(Nb,N1,N1,3) :: xtc       !panel chordwise tangents
real,    intent(in), dimension(Nb,N1,N1,3) :: xts       !panel spanwise tangents
complex, intent(in), dimension(Nb,N1,N1,3) :: dxc       !collocation points
complex, intent(in), dimension(Nb,N1,N1,3) :: dxn       !panel normals

!outputs
complex, intent(out), dimension(N2,N2) :: da          !influence matrix
complex, intent(out), dimension(N2)   :: drhs        !rhs of linear system

real, intent(out), dimension(N2,N2) :: daGc, daGs   !influence matrix for tangent velocities due to every other vortex ring

!internal variables
integer :: i,j, p,q, r,s, t
integer :: countC,countV
real, dimension(3) :: u
complex :: v
complex :: phase0, phase, shift

!initialise
u    = 0.0
v    = 0.0
da   = 0.0
drhs = 0.0

phase0 = -imag*omega*dt
!phase0 = 1.0-exp(-imag*omega*dt)
shift  = exp(imag*omega*dt)

!influence matrix arranged in blocks of the blade section at each spanwise location
!scan over collocation points
do 1 s = 1,Nb
do 1 j = 1,Ns
do 1 i = 1,Nc
        
        countC = indx3(i,j,Nc,s,Ns)
        
        !u.d(xn), and (dxc/dt).xn
        drhs(countC) = -imag*omega*dot_product(dxc(s,i,j,:) ,  xn(s,i,j,:)) &
                                 - dot_product(uinf         , dxn(s,i,j,:))
        
        do 1 r = 1,Nb
        do 1 q = 1,Ns

                !scan over bound vortex rings
                do 2 p = 1,Nc
                        
                        countV = indx3(p,q,Nc,r,Ns)
                        
                        u = vring(xc(s,i,j,:), xv(r,p  ,q  ,:), xv(r,p  ,q+1,:), &
                                               xv(r,p+1,q+1,:), xv(r,p+1,q  ,:))
                        
                        da(countC,countV) = da(countC,countV) + dot_product(u , xn(s,i,j,:))
                        
                        daGc(countC,countV) = daGc(countC,countV) + dot_product(u, xtc(s,i,j,:))
                        daGs(countC,countV) = daGs(countC,countV) + dot_product(u, xts(s,i,j,:))
                        
                2 continue
                
                !scan over shed vortex rings
                phase = phase0
                do 1 p = Nc+1,Nw-1
                        
                        u = vring(xc(s,i,j,:), xv(r,p  ,q  ,:) , xv(r,p  ,q+1,:), &
                                               xv(r,p+1,q+1,:) , xv(r,p+1,q  ,:))
                        
                        v = phase*dot_product(u , xn(s,i,j,:))
                        
                        do 3 t = 1,Nc
                                
                                countV = indx3(t,q,Nc,r,Ns)
                                
                                da(countC,countV) = da(countC,countV) + v
                                
                                daGc(countC,countV) = daGc(countC,countV) + dot_product(u, xtc(s,i,j,:))
                                daGs(countC,countV) = daGs(countC,countV) + dot_product(u, xts(s,i,j,:))
                                
                        3 continue
                        
                        phase = phase*shift
                        
1 continue

do 4 s = 1,Nb
do 4 j = 1,Ns
do 4 i = 1,Nc
        
        countC = indx3(i,j,Nc,s,Ns)
        
        u = vring(xc(s,i,j,:), xv(s,i  ,j  ,:), xv(s,i  ,j+1,:), &
                               xv(s,i+1,j+1,:), xv(s,i+1,j  ,:))
        
        daGc(countC,countC) = daGc(countC,countC) - dot_product(u, xtc(s,i,j,:))
        daGs(countC,countC) = daGs(countC,countC) - dot_product(u, xts(s,i,j,:))

4 continue
        
 
end subroutine flex_matrix
