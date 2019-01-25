subroutine nbody_deallocate_arrays
! Code deallocates nbody arrays ready for the next disc model

use embryodata
implicit none

deallocate(pos,vel,acc)
deallocate(newpos,newvel)
deallocate(mass,angmom,angmag)
deallocate(ekin,epot,etot)


end subroutine nbody_deallocate_arrays
