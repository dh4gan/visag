subroutine nbody_acceleration(position,velocity,acceleration)
! Routine drives calculation of all different acceleration terms

use embryodata,only: nbodies
implicit none


real,dimension(3,nbodies), intent(in) :: position,velocity
real,dimension(3,nbodies), intent(out) :: acceleration

acceleration(:,:) = 0.0

call nbody_grav_acceleration(position,acceleration)
call nbody_drag_terms(position,velocity,acceleration)

end subroutine nbody_acceleration
