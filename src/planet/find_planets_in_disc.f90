subroutine find_planets_in_disc
!
! Finds grid cell planets are in
!

use gravdata, only: rz,isr,nrgrid
use planetdata
implicit none

integer :: icheck,iplanet

icheck = 0
do iplanet=1,nplanet

      ! Check from vicinity of last known location
      iplanetrad(iplanet) = iplanetrad(iplanet)-int(0.1*nrgrid)

      ! If at start or near inner boundary, use it as starting point
      if(iplanetrad(iplanet)<isr) iplanetrad(iplanet) = isr

   do while(rz(iplanetrad(iplanet))<ap(iplanet))
      iplanetrad(iplanet) = iplanetrad(iplanet)+1
   enddo

enddo


end subroutine find_planets_in_disc
