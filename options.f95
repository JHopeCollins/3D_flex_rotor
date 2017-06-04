!set parameters for flexible 3d blade

module options

use constants

!variables
implicit none

!options

!programme sections
integer, parameter :: Seigen = 1
integer, parameter :: Sloads = 1
integer, parameter :: Ueigen = 0
integer, parameter :: Uloads = 0

!blade geometry
real, parameter :: rootchord =  1.0       !chord length at blade root
real, parameter :: tipchord  =  0.4       !chord length at blade tip
real, parameter :: leadaxis  =  0.2       !relative size of streamwise axis of le ellipsoid
real, parameter :: AR        =  6.0       !aspect ratio relative to root chord
real, parameter :: clip      =  1.00      !clipped wing tips for ellipsoid

!rotor geometry
integer, parameter :: Nb     = 3          !number of blades
real,    parameter :: pitch  = 30.0*d2r   !pitch angle of blades
real,    parameter :: sweep  =-00.0*d2r   !sweep angle (+ve for leaning streamwise)
real,    parameter :: TSR    = 4.0        !tip speed ratio
real,    parameter :: hub    = 1.0        !hub radius

!camberline geometry
real, parameter :: h = 0.30             !maximum camber height
real, parameter :: m = 0.40             !maximum camber location

!freestream velocity
!real, parameter :: alpha = 0.0*d2r
!real, parameter, dimension(3) :: uinf = [cos(alpha) , 0.0 , sin(alpha)]
real, parameter, dimension(3) :: uinf = [1.0 , 0.0 , 0.0]

!flexing wavenumbers
real, parameter :: lamda = 3.0                  !wavelength relative to chord
real, parameter :: k     = pi/lamda             !chordwise reduced frequency
real, parameter :: omega = 2.0*k                !frequency / wavenumber

!discretisation values
integer, parameter :: Nwp = 10                          !number of wake wavelengths
integer, parameter :: res = 10                          !wake panels/wavelength
real,    parameter :: panelAR = 5.0                     !panel aspect ratio

integer, parameter :: Nc = int(res/lamda)               !number of chordwise panels
integer, parameter :: Ns = Nc*ceiling(AR/panelAR)       !number of spanwise panels
integer, parameter :: Nw = Nc + Nwp*res                 !number of last wake panel
real,    parameter :: dt  = lamda/res

!derived values
real, parameter :: span = rootchord*AR           !blade spanwise length
real, parameter :: trailaxis =  1.0-leadaxis     !relative length of streamwise axis for te ellipsoid
real, parameter :: Romega = TSR/(hub+span)       !rotor angular velocity

!blade collocation points in local frame of reference
real, dimension(N1,N1,3) :: blade

end module  options
