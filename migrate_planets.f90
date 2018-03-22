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

! Assume that planet torque contributes to its migration only
! (Detailed torque balance between disc and individual bodies)


do iplanet=1,nplanet

! Integrate torque*sigma over the entire disc

    do i=isr,ier
    adot(iplanet) = adot(iplanet) + torquei(iplanet,i)*sigma(i)*(rz(i)-rz(i-1))
    enddo

    ! Multiply by appropriate factors to get adot

    adot(iplanet) = -adot(iplanet)*4.0*pi*ap(iplanet)/mp(iplanet)

    ! get overall migration timescale
    tmig(iplanet) = ap(iplanet)/abs(adot(iplanet))

    ap(iplanet) = ap(iplanet) + adot(iplanet)*dt

enddo


end subroutine migrate_planets