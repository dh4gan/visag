subroutine setup_planets(planetfile)
!
! Subroutine sets up planets to be added to disc
! Planet data read from separate input file
!
!

use planetdata
character(100), intent(in) :: planetfile

! Open planet file

open(10, file=planetfile,status='old')
! Read input data

read(10,*) nplanet

nactive = nplanet

allocate(mp(nplanet),ap(nplanet))
allocate(lambdaI(nplanet,nmax), lambdaII(nplanet,nmax))
allocate(fII(nplanet,nmax))
allocate(adot(nplanet),tmig(nplanet),tmigI(nplanet))

mp(:) = 0.0
ap(:) = 0.0
lambdaII(:,:) = 0.0
lambdaI(:,:) = 0.0
fII(:,:) = 0.0
adot(:) = 0.0
tmig(:) = 0.0
tmigI(:) = 0.0

do iplanet=1,nplanet
read(10,*) mp(iplanet), ap(iplanet)
enddo



end subroutine setup_planets