subroutine nbody_rk4
! This subroutine drives the N-Body integration over a single timestep
! Integration is done in a separate set of arrays for speed
! Orbital elements are stored in GE_embryo type

! Do integration

use stardata,only: debug
use embryodata

implicit none

integer :: ibody
logical :: withintolerance

withintolerance = .false.
if(debug=='y') print*, 'Attempting integration RK4'

do while(withintolerance .eqv. .false.)

newpos(:,:) = 0.0
newvel(:,:) = 0.0

! Switch off velocities of finished particles
do ibody=2,nbodies
   if(embryo(ibody-1)%finished==1) vel(:,ibody)=0.0
enddo

call nbody_integrate(dt_nbody,pos,vel,newpos,newvel)
call nbody_timestep(newpos,newvel)

if(maxerror<tolerance) withintolerance=.true.

end do

pos = newpos
vel = newvel

call nbody_system_properties

end subroutine nbody_rk4
