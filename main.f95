!panel method model for 3d thin flexible blade with small deformation and NACA 2 digit camberline and
!prescribed wake

program main

use constants
use calls
use options

!variables
implicit none

real, dimension(Nb,N1,N1,3) :: x             !grid points (panel corners)
real, dimension(Nb,N1,N1,3) :: xc            !panel collocation points
real, dimension(Nb,N1,N1,3) :: xv            !vortex ring corners
real, dimension(Nb,N1,N1,3) :: xn            !panel normals
real, dimension(Nb,N1,N1,3) :: xtc           !panel tangents chordwise
real, dimension(Nb,N1,N1,3) :: xts           !panel tangents spanwise

complex, dimension(Nb,N1,N1,3) :: dx         !grid points (panel corners)
complex, dimension(Nb,N1,N1,3) :: dxc        !panel collocation points
complex, dimension(Nb,N1,N1,3) :: dxv        !vortex ring corners
complex, dimension(Nb,N1,N1,3) :: dxn        !panel normals

real, dimension(N2,N2) :: a, adum            !influence matrix
real, dimension(N2)   :: rhs                 !free stream influence

real, dimension(N2,N2) :: aGc, aGs           !influence matrices for chordwise and spanwise tangential velocities (global)
real, dimension(N2,N2) :: alc, als           !influence matrices for chordwise and spanwise tangential velocities (local)
real, dimension(N2,N2) :: daGc, daGs         !influence matrices for chordwise and spanwise tangential velocities
real, dimension(Nb,N1,N1) :: SGc, SGs        !chordwise and spanwise velocity due to every other panel (global)
real, dimension(Nb,N1,N1) :: Slc, Sls        !chordwise and spanwise velocity self induced by panel (local)

real, dimension(Nb,N1,N1) :: UGc, UGs        !chordwise and spanwise velocity due to every other panel (global)
real, dimension(Nb,N1,N1) :: Ulc, Uls        !chordwise and spanwise velocity self induced by panel (local)

complex, dimension(N2,N2) :: da, dadum       !influence matrix
complex, dimension(N2)    :: drhs            !free stream influence
complex, dimension(N2,N2) :: dldg            !unsteady lift per circulation

real   , dimension(Nb,N1,N1) :: g            !bound vorticity distribution
complex, dimension(Nb,N1,N1) :: dg           !bound vorticity distribution
real :: lift                                 !lift on blade
real :: drag                                 !wake induced drag
real :: torque

integer :: i, j, s, p

!lapack variables
integer :: info1, info2, lda = (Nb*Ns*Nc)
integer, dimension(N2)   :: ipiv
complex, dimension(5*N2) :: CWORK         !cgeev workspace
real, dimension(2*N2)    :: RWORK         !cgeev workspace

!eigen variables
complex, dimension(N2)   :: deigval       !eigenvalues of da
complex, dimension(N2,N2) :: deigvecL     !left  eigenvector i of da is deigvecL(:,i)
complex, dimension(N2,N2) :: deigvecR     !right eigenvector i of da is deigvecR(:,i)

real, dimension(5*N2)    :: SWORK         !cgeev workspace
real, dimension(N2)    :: eigval_real     !real part of eigenvalues of a
real, dimension(N2)    :: eigval_imag     !imag part of eigenvalues of a
complex, dimension(N2) :: eigval          !eigenvalues of a
real, dimension(N2,N2)  :: eigvecL        !left  eigenvector i of a is eigvecL(:,i)
real, dimension(N2,N2)  :: eigvecR        !right eigenvector i of a is eigvecR(:,i)

open(unit=4, file='Slefteigvec.dat', action='write', status='replace')
open(unit=3, file='Sgam.dat',        action='write', status='replace')
open(unit=2, file='gam.dat',         action='write', status='replace')

!initialise
        i = 1
        j = 1
        p = 1
        s = 1
        a = 0.0
        g = 0.0
        x = 0.0
        xn = 0.0
        xc = 0.0
        xv = 0.0
        da = 0.0
        dg = 0.0
        dx = 0.0
        dxn = 0.0
        dxc = 0.0
        dxv = 0.0
        rhs = 0.0
        drhs = 0.0
        lift = 0.0
        drag = 0.0
        dadum = 0.0
        torque = 0.0

print *, 'Nb', Nb, 'Ns', Ns, 'Nc', Nc, 'Nw', Nw-Nc, 'lda', lda
!---------steady solution-------------------------------------------------------


!panel grid
call grid(x,xc,xv, xn,xtc,xts)
print *, 'grid ok'

if ((Seigen.eq.1) .or. (Sloads.eq.1) .or. (Ueigen.eq.1) .or. (Uloads.eq.1)) then
!influence matrix and rhs
call matrix(x,xc,xv, xn,xtc,xts, a,rhs, aGc,aGs, alc,als)
print *, 'matrix ok'
end if


!---------eigenvectors--------------
if (Seigen.eq.1) then
adum = a

call sgeev('V', 'V', lda, a, N2, eigval_real, eigval_imag, eigvecL, N2, eigvecR, N2, SWORK, 5*N2, info1)

print *, 'steady eig info', info1
eigval = cmplx(eigval_real, eigval_imag)

do 96 p = 1,10
        do 93 j = 1,Ns
        do 93 i = 1,Nc
        !write(1000+p,*) blade(i,j,1:2), eigvecL(indx3(i,j,Nc,1,Ns),lda+1-p)
        write(1000+p,*) blade(i,j,1:2), eigvecL(indx3(i,j,Nc,1,Ns),p)
        93 continue
96 continue

a = adum
end if
!---------eigenvectors--------------


!---------circulations--------------
if (Sloads.eq.1) then

call sgetrf(     lda, lda, a, N2, ipiv,         info1)
call sgetrs('n', lda, 1,   a, N2, ipiv, rhs, N2, info2)

print *, 'steady info1', info1
print *, 'steady info2', info2

do 90 s = 1,Nb
        do 90 j = 1,Ns
        do 90 i = 1,Nc
        g(s,i,j) = rhs(indx3(i,j,Nc,s,Ns))
90 continue

do 91 j = 1,Ns
        write(3,*) blade(1,j,1:2), (g(1,1,j)), g(1,1,j)
        do 92 i = 2,Nc
        write(3,*) blade(i,j,1:2), g(1,i,j), g(1,i,j)-g(1,i-1,j)
        92 continue
91 continue

Slc = rmatprod2mat(alc,g)
Sls = rmatprod2mat(als,g)
SGc = rmatprod2mat(aGc,g)
SGs = rmatprod2mat(aGs,g)

!call loads(x, g, lift, drag)

end if
!---------circulations--------------


!-------unsteady solution-------------------------------------------------------
if ((Ueigen.eq.1) .or. (Uloads.eq.1)) then

!change in panel grid
call flex_grid(x,xn, dx,dxc,dxv, dxn)
print *, 'flex grid ok'


!influence matrix and rhs for flexible oscillations
call flex_matrix(xc,xv, xn,xtc,xtc, dxc,dxn, da,drhs, daGc,daGs)
print *, 'flex matrix ok'

!unsteady lift per circulation matrix
call flex_liftmatrix(x,xc,xv, xtc,xts, SGc,SGs, Slc,Sls, dldg)
print *, 'flex lift matrix ok'

end if


!---------eigenvectors--------------
if (Ueigen.eq.1) then
dgdu = da

call cgeev('V', 'V', lda, da, N2, deigval, deigvecL, N2, deigvecR, N2, CWORK, 5*N2, RWORK, info1)

print *, 'flex eig info', info1

do 99 p = 1,10
        do 99 j = 1,Ns
        do 99 i = 1,Nc
        write(2000+p,*) blade(i,j,1:2), real(deigvecL(indx3(i,j,Nc,s,Ns),p)), aimag(deigvecL(indx3(i,j,Nc,s,Ns),p))
99 continue

da = dgdu
end if
!---------eigenvectors--------------


!---------circulations--------------
if (Uloads.eq.1) then

call cgetrf(     lda, lda, da, N2, ipiv,              info1)
call cgetri(     lda,      da, N2, ipiv, CWORK, 5*N2, info2)
!call cgetrs('n', lda, 1,   da, N2, ipiv, drhs, N2,    info2)

print *, 'flexible info1', info1
print *, 'flexible info2', info2
dgdu = da

do 80 s = 1,Nb
        do 80 j = 1,Ns
        do 80 i = 1,Nc
        dg(s,i,j) = drhs(indx3(i,j,Nc,s,Ns))
80 continue


do 97 j = 1,Ns
        write(2,*) xc(1,1,j,1), xc(1,1,j,3), real(dg(1,1,j)), aimag(dg(1,1,j)), real(dg(1,1,j)), aimag(dg(1,1,j))
        do 98 i = 2,Nc
        write(2,*) xc(1,i,j,1), xc(1,i,j,3), real(dg(1,i,j)), aimag(dg(1,i,j)) &
                , real(dg(1,i,j)-dg(1,i-1,j)), aimag(dg(1,i,j)-dg(1,i-1,j))
        98 continue
97 continue

end if
!---------circulations--------------



!-------functions---------------------------------------------------------------
contains

function rmatprod2mat(a,b)
        use constants
        use options
        use calls
        
        !variables
        implicit none
        
        !inputs
        real, dimension(N2,N2) :: a
        real, dimension(N2)   :: b
        
        !outputs
        real, dimension(Nb,N1,N1) :: rmatprod2mat
        
        !internal variables
        !real, dimension(N2) :: vec
        integer :: i, j, s
        
        rmatprod2mat = 0.0
        
        b = matmul(a,b)
        
        do 1 s = 1,Nb
        do 1 j = 1,Ns
        do 1 i = 1,Nc
                rmatprod2mat(s,i,j) = b(indx3(i,j,Nc,s,Ns))
        1 continue
end function rmatprod2mat


function cmatprod2mat(a,b)
        use constants
        use options
        use calls
        
        !variables
        implicit none
        
        !inputs
        complex, dimension(N2,N2) :: a
        complex, dimension(N2)   :: b
        
        !outputs
        complex, dimension(Nb,N1,N1) :: cmatprod2mat
        
        !internal variables
        !complex, dimension(N2) :: vec
        integer :: i, j, s
        
        cmatprod2mat = 0.0
        
        b = matmul(a,b)
        
        do 1 s = 1,Nb
        do 1 j = 1,Ns
        do 1 i = 1,Nc
                cmatprod2mat(s,i,j) = b(indx3(i,j,Nc,s,Ns))
        1 continue
end function cmatprod2mat


end program main
