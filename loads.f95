!load calculation for 3d thin wing under steady free stream

subroutine loads(x, g, lift, drag)

use constants
use calls
use options

!variables
implicit none

!inputs
real, intent(in), dimension(Nb,N1,N1,3) :: x
real, dimension(Nb,N1,N1)   :: g

!outputs
real, intent(out) :: lift
real, intent(out) :: drag

!internal variables
real, dimension(Nb,N1,N1) :: dl, dd
real :: dy
integer :: i, j, s

!initialise
lift = 0.0
drag = 0.0
dl = 0.0
dd = 0.0
dy = 0.0

do 1 s = 1,Nb
do 1 j = 1,Ns
        
        dy =      x(s,1,j+1,2) - x(s,1,j,2)
        dy = dy + x(s,2,j+1,2) - x(s,2,j,2)
        dy = dy /2.0
        
        dl(Nb,1,j) = g(Nb,1,j)*dy

1 continue

do 2 s = 1,Nb
do 2 j = 1,Ns
do 2 i = 2,Nc
        
        dy =      x(s,i  ,j+1,2) - x(s,i  ,j,2)
        dy = dy + x(s,i+1,j+1,2) - x(s,i+1,j,2)
        dy = dy /2.0
        
        dl(s,i,j) = (g(s,i,j)- g(s,i-1,j))*dy

2 continue

lift = sum(dl(1:Nb,1:Nc,1:Ns))
drag = sum(dd(1:Nb,1:Nc,1:Ns))


end subroutine loads


!calculate the unsteady lift per unit unsteady circulation matrix

subroutine flex_liftmatrix(x,xc,xv, xtc,xts, SGc,SGs, Slc,Sls, dldg)

use constants
use calls
use options

!variables
implicit none

!inputs
real, intent(in), dimension(Nb,N1,N1,3) :: x
real, intent(in), dimension(Nb,N1,N1,3) :: xc
real, intent(in), dimension(Nb,N1,N1,3) :: xv
real, intent(in), dimension(Nb,N1,N1,3) :: xtc
real, intent(in), dimension(Nb,N1,N1,3) :: xts
real, intent(in), dimension(Nb,N1,N1)   :: SGc
real, intent(in), dimension(Nb,N1,N1)   :: SGs
real, intent(in), dimension(Nb,N1,N1)   :: Slc
real, intent(in), dimension(Nb,N1,N1)   :: Sls

!outputs
complex, intent(out), dimension(N2,N2) :: dldg

!internal variables
integer :: i,j, p,q, r,s, t
integer :: countC, countV
real, dimension(3) :: u
real :: dx, dy, ds
complex :: phase0, phase, shift

!initialise
u = 0
dldg = 0.0
phase0 = -imag*omega*dt
phase0 = 1.0-exp(-imag*omega*dt)
shift  = exp(imag*omega*dt)


do 1 s = 1,Nb
do 1 j = 1,Ns
do 1 i = 1,Nc
        
        countC = indx3(i,j,Nc,s,Ns)
        
        do 1 r = 1,Nb
        do 1 q = 1,Ns
                
                !bound vortices
                do 2 p = 1,Nc
                        
                        countV = indx3(p,q,Nc,r,Ns)
                        
                        u = vring(xc(s,i,j,:), xv(r,p  ,q  ,:), xv(r,p  ,q+1,:), &
                                               xv(r,p+1,q+1,:), xv(r,p+1,q  ,:))
                        
                        dldg(countC,countV) = dldg(countC,countV) + Slc(s,i,j)*dot_product(u,xtc(s,i,j,:))
                        dldg(countC,countV) = dldg(countC,countV) + Sls(s,i,j)*dot_product(u,xts(s,i,j,:))
                        
                2 continue
                
                !shed vortex rings
                phase = phase0
                do 1 p = Nc+1,Nw-1
                        
                        u = vring(xc(s,i,j,:), xv(r,p  ,q  ,:), xv(r,p  ,q+1,:), &
                                               xv(r,p+1,q+1,:), xv(r,p+1,q  ,:))
                        
                        do 3 t = 1,Nc
                                
                                countV = indx3(t,q,Nc,r,Ns)
                                
                                dldg(countC,countV) = dldg(countC,countV) + phase*Slc(s,i,j)*dot_product(u,xtc(s,i,j,:))
                                dldg(countC,countV) = dldg(countC,countV) + phase*Sls(s,i,j)*dot_product(u,xts(s,i,j,:))
                        
                        3 continue
                        
                        phase = phase*shift
                        
1 continue

do 4 s = 1,Nb
do 4 j = 1,Ns
do 4 i = 1,Nc
        
        countC = indx3(i,j,Nc,s,Ns)
        
        dy = norm2(0.5*(x(s,i  ,j+1,:) - x(s,i  ,j,:) + &
                        x(s,i+1,j+1,:) - x(s,i+1,j,:)))
        
        dx = norm2(0.5*(x(s,i+1,j  ,:) - x(s,i,j  ,:) + &
                        x(s,i+1,j+1,:) - x(s,i,j+1,:)))
        
        u = vring(xc(s,i,j,:), xv(s,i  ,j  ,:), xv(s,i  ,j+1,:), &
                               xv(s,i+1,j+1,:), xv(s,i+1,j  ,:))
        
        dldg(countC,countC) = dldg(countC,countC) - Slc(s,i,j)*dot_product(u, xtc(s,i,j,:))
        dldg(countC,countC) = dldg(countC,countC) - Sls(s,i,j)*dot_product(u, xts(s,i,j,:))
        
        dldg(countC,countC) = dldg(countC,countC) + SGc(s,i,j)/dx
        dldg(countC,countC) = dldg(countC,countC) + SGs(s,i,j)/dy
        
        if (i.ne.1) then
                
                dldg(countC,countC-1)  = dldg(countC,countC-1)  - SGc(s,i,j)/dx
        end if

        if (j.ne.1) then
                
                dldg(countC,countC-Nc) = dldg(countC,countC-Nc) - SGs(s,i,j)/dy
        end if
        
4 continue

do 5 s = 1,Nb
        do 5 j = 1,Ns
        do 5 i = 1,Nc
                
        countC = indx3(i,j,Nc,s,Ns)

        ds = norm2(rcross(x(s,i+1,j  ,:)-x(s,i,j,:), &
                          x(s,i  ,j+1,:)-x(s,i,j,:)))
        
        dldg(countC,:) = dldg(countC,:)*ds
        
5 continue


end subroutine flex_liftmatrix
