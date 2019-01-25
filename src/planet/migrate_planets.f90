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

!
! Assume that planet torque contributes to its migration only
! (Detailed torque balance between disc and individual bodies)
!

do iplanet=1,nplanet

   ! Integrate torque*sigma over the entire disc
   adot(iplanet) = 0.0
    
   do i=isr,ier       
      adot(iplanet) = adot(iplanet) + torquei(iplanet,i)*sigma(i)/drzm1(i)
   enddo

    ! Multiply by appropriate factors to get adot    
   adot(iplanet) = -adot(iplanet)*4.0*pi* &
        omegaK(iplanetrad(iplanet))*ap(iplanet)*ap(iplanet)/mp(iplanet)
   !adot(iplanet) = -adot(iplanet)*4.0*pi*ap(iplanet)/mp(iplanet)
    
    ! get numerically determined migration timescale

    if(abs(adot(iplanet))>1.0e-30) then
       tmig(iplanet) = ap(iplanet)/(-adot(iplanet))
    else
       tmig(iplanet) = 0.0
    endif


    ! TODO - write new routine to choose between basic and N-Body motion
    
    ! Move planets
    ap(iplanet) = ap(iplanet) + adot(iplanet)*dt
  
enddo

end subroutine migrate_planets
