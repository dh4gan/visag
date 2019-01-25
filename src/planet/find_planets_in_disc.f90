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

      do while(rz(iplanetrad(iplanet))<ap(iplanet) .and. iplanetrad(iplanet)<nrgrid)
         iplanetrad(iplanet) = iplanetrad(iplanet)+1
      enddo

      ! If planet beyond outer disc radius, then set this to a negative value
      ! Allows us to skip torque calculations later
      
      if(iplanetrad(iplanet)>nrgrid) iplanetrad(iplanet) = -10

      ! Planets that move inside the minimum radius are removed here
      if(rz(iplanetrad(iplanet)) < rremove) alive(iplanet)=0
      
   
enddo


end subroutine find_planets_in_disc
