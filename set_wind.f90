SUBROUTINE set_wind
! Routine takes inputs and sets up disc winds

use gravdata
use winddata
use unitdata

  ! Convert wind mass loss rate into a surface density loss rate
  ! Assume wind operates from rwind -> 25 AU for normalization
  ! Use sigmadot proportional to 1/r 

  rwindout = 25.0d0*AU

  if ((rwindout-rwind) .gt. 0.0d0) then 
     sigdot = mdot_wind / (2.0d0*pi*(rwindout-rwind))
  else
     sigdot = 0.0d0
  endif

return
END SUBROUTINE set_wind
