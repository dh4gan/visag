subroutine migrate_planets
!
! Subroutine uses the torque exerted by each planet to find the
! migration timescale for each planet
!

use planetdata
use gravdata
use unitdata

implicit none

integer :: iplanet,i
real :: tmigcheck,tmigcheck2

! Assume that planet torque contributes to its migration only
! (Detailed torque balance between disc and individual bodies)


do iplanet=1,nplanet

! Integrate torque*sigma over the entire disc
   adot(iplanet) = 0.0
    
   do i=isr,ier       
      adot(iplanet) = adot(iplanet) + torquei(iplanet,i)*sigma(i)/drzm1(i)
   enddo

    ! Multiply by appropriate factors to get adot

    !adot(iplanet) = -adot(iplanet)*4.0*pi*G*mstar/(omegaK(iplanetrad(iplanet))*mp(iplanet)*ap(iplanet))
   adot(iplanet) = -adot(iplanet)*4.0*pi*G/(omegaK(iplanetrad(iplanet))*ap(iplanet))
    
    call calc_typeI_migration(iplanet,tmigcheck)
    ! get overall migration timescale
    tmig(iplanet) = ap(iplanet)/abs(adot(iplanet))
       
    !print*, ap(iplanet)/AU, tmig(iplanet)/yr, tmigcheck/yr

    if(tmigcheck<0.0) then
       print*, 'NEGATIVE TMIG1'
       stop
       endif

    ap(iplanet) = ap(iplanet) + adot(iplanet)*dt

enddo

!STOP
end subroutine migrate_planets
