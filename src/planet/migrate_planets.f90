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
real :: tmigcheck,tmigcheck2,hp,aspectratio,mdiscmig

! Assume that planet torque contributes to its migration only
! (Detailed torque balance between disc and individual bodies)


do iplanet=1,nplanet

! Integrate torque*sigma over the entire disc
   adot(iplanet) = 0.0
    do i=isr,ier
       
    adot(iplanet) = adot(iplanet) + torquei(iplanet,i)*sigma(i)/drzm1(i)
    enddo

    ! Multiply by appropriate factors to get adot

    !adot(iplanet) = -adot(iplanet)*4.0*pi*ap(iplanet)/mp(iplanet)
    adot(iplanet) = -adot(iplanet)*4.0*pi*G*sqrt(ap(iplanet)/(G*Mstar))

    ! get overall migration timescale
    tmig(iplanet) = ap(iplanet)/abs(adot(iplanet))

    hp = cs(iplanetrad(iplanet))/omegaK(iplanetrad(iplanet))
    aspectratio = hp/ap(iplanet)

    mdiscmig = pi*ap(iplanet)*ap(iplanet)*sigma(iplanetrad(iplanet))

    tmigcheck = mstar*mstar*aspectratio*aspectratio/(mdiscmig*mp(iplanet)*omegaK(iplanetrad(iplanet)))

    tmigcheck2 = mstar*cs(iplanetrad(iplanet))/(mp(iplanet)*ap(iplanet)*omegaK(iplanetrad(iplanet))**2)

    
   print*, iplanet, adot(iplanet), ap(iplanet), tmig(iplanet)/yr, tmigcheck/yr, tmigcheck2/yr      
        

    ap(iplanet) = ap(iplanet) + adot(iplanet)*dt

enddo

STOP

end subroutine migrate_planets
