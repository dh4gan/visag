subroutine move_planets
  !
  ! Subroutine uses the torque exerted by each planet to find the
  ! migration timescale for each planet
  !
  ! If N Body integration is requested, N body integrator called here
  !

use planetdata
use gravdata
use unitdata
use nbodydata, only: a

implicit none

integer :: iplanet,i

!
! Assume that planet torque contributes to its migration only
! (Detailed torque balance between disc and individual bodies)
!

do iplanet=1,nplanet

   ! Integrate torque*sigma over the entire disc
   adot(iplanet) = 0.0

   if(iplanetrad(iplanet)>0) then
      do i=isr,ier       
         adot(iplanet) = adot(iplanet) + torquei(iplanet,i)*sigma(i)/drzm1(i)
      enddo

      ! Multiply by appropriate factors to get adot    
      adot(iplanet) = -adot(iplanet)*4.0*pi* &
           omegaK(iplanetrad(iplanet))*ap(iplanet)*ap(iplanet)/mp(iplanet)
      !adot(iplanet) = -adot(iplanet)*4.0*pi*ap(iplanet)/mp(iplanet)
   endif
   
    ! get numerically determined migration timescale

    if(abs(adot(iplanet))>1.0e-30) then
       tmig(iplanet) = ap(iplanet)/(-adot(iplanet))
    else
       tmig(iplanet) = 0.0
    endif

 enddo

 ! Migration timescales calculated - now must move bodies

 ! If using N Body integrator, call it here
 
 if(nbodychoice=="y") then

    call nbody_rk4

    ! Update planet data after integration step
    do iplanet=1,nplanet
       ap(iplanet) = a(iplanet+1)/AU
    enddo
    
 else

    ! If not using an integrator, simply move planets radially
    
    do iplanet=1,nplanet   
       ap(iplanet) = ap(iplanet) + adot(iplanet)*dt
    enddo
       
 endif
    

end subroutine move_planets
