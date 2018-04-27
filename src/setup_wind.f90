SUBROUTINE setup_wind
!
! Routine takes inputs and sets up disc winds
! X Ray winds from Owen et al (2010)
! FUV/EUV winds from Pont et al, Alexander et al
!

use gravdata
use winddata
use unitdata

implicit none

integer :: i
real :: xwind,sigdotxray


allocate(sigdot_wind(nmax))

sigdot_wind(:) = 0.0

print*, 'Setting up wind parameters'

mdot_wind = 0.0

do i = isr, ier
xwind = 0.85*(rz(i)/1.496d13)*(mstar/1.989d33)**(-1.0d0)

call xraywind(sigdotxray,xwind)

mdot_wind = mdot_wind + sigdotxray*2.0d0*3.14159*rz(i)/drzm1(i)

enddo

windnorm = 6.25d-9*(mstar/solarmass)**(-0.068)*(Lx/1.0d30)**1.14/mdot_wind

return
END SUBROUTINE setup_wind
