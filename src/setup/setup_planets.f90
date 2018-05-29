subroutine setup_planets
!
! Subroutine sets up planets to be added to disc
! Planet data read from separate input file
!
!

use gravdata
use planetdata
use unitdata

implicit none

integer :: iplanet

! Open planet file

open(10, file=planetfile,status='old')
! Read input data

read(10,*) nplanet

print*, 'There are ',nplanet, 'planets'
nactive = nplanet

allocate(mp(nplanet),ap(nplanet),alive(nplanet), iplanetrad(nplanet))
allocate(lambdaI(nplanet,nmax), lambdaII(nplanet,nmax))
allocate(fII(nplanet,nmax))
allocate(adot(nplanet),tmig(nplanet),tmigI(nplanet))
allocate(torquei(nplanet,nmax), torque_term(nmax), total_planet_torque(nmax))


alive(:) = 1
mp(:) = 0.0
ap(:) = 0.0
iplanetrad(:) = 0

lambdaII(:,:) = 0.0
lambdaI(:,:) = 0.0
fII(:,:) = 0.0
adot(:) = 0.0
tmig(:) = 0.0
tmigI(:) = 0.0
torquei(:,:) = 0.0
total_planet_torque(:) = 0.0

do iplanet=1,nplanet
   read(10,*) mp(iplanet), ap(iplanet)
enddo
 

! Convert to correct units
mp(:) = mp(:)*mjup
ap(:) = ap(:)*AU

do iplanet=1,nplanet
  ! Find the location of each planet in the disc
   iplanetrad(iplanet) = isr

   do while (rz(iplanetrad(iplanet))<ap(iplanet))
      iplanetrad(iplanet) = iplanetrad(iplanet)+1
   enddo

   print*, 'Planet ', iplanet, ' located at cell ', iplanetrad(iplanet)
   print*, 'Radius: ', rz(iplanetrad(iplanet))/AU

enddo

end subroutine setup_planets
