subroutine nbody_system_properties
! Calculates the various properties of the system
! Includes energy, angular momentum and orbits

use embryodata

implicit none

call angular_momentum
call energy

totalmass = sum(mass)
call calc_orbit_from_vector(totalmass)

end subroutine nbody_system_properties

subroutine angular_momentum
! Calculate the angular momentum of all the bodies

use embryodata

implicit none

integer :: ix

! Compute angular momentum = r x v

angmom(1,:) = pos(2,:)*vel(3,:) - pos(3,:)*vel(2,:)
angmom(2,:) = pos(3,:)*vel(1,:) - pos(1,:)*vel(3,:)
angmom(3,:) = pos(1,:)*vel(2,:) - pos(2,:)*vel(1,:)
angmag(:) = sqrt(angmom(1,:)*angmom(1,:) + angmom(2,:)*angmom(2,:) + angmom(3,:)*angmom(3,:))

do ix=1,3
system_angmom(ix) = sum(angmom(ix,:))
enddo

system_ang =sqrt(system_angmom(1)*system_angmom(1) + system_angmom(2)*system_angmom(2)+system_angmom(3)*system_angmom(3))

if(initial_system_ang>1.0e-30)then
dL = (system_ang-initial_system_ang)/initial_system_ang
else
dL=0.0
endif

end subroutine angular_momentum


subroutine energy
! Calculate the orbital parameters of all the bodies
! Also compute energy and angular momentum for error tracking

use embryodata

implicit none

integer :: ix, ibody,jbody
real :: relpos
real,dimension(nbodies) :: vmag
real,dimension(3) :: sep

! Compute kinetic energy of the bodies

vmag(:) = sqrt(vel(1,:)*vel(1,:) + vel(2,:)*vel(2,:) + vel(3,:)*vel(3,:))
ekin(:) = 0.5*vmag(:)*vmag(:)

! Compute potential energy

do ibody=1,nbodies

epot(ibody)=0.0

do jbody=1,nbodies

if(ibody==jbody) cycle

do ix=1,3
sep(ix) = pos(ix,ibody) - pos(ix,jbody)
enddo

relpos = sqrt(sep(1)*sep(1) + sep(2)*sep(2)+sep(3)*sep(3) +rsoft*rsoft)

do ix=1,3
epot(ibody) = epot(ibody) - mass(jbody)/(relpos)
enddo
enddo

enddo

etot(:) = ekin(:) + epot(:)

system_energy = sum(etot)

if(initial_system_energy >1.0e-30) then
dE = (system_energy-initial_system_energy)/initial_system_energy
else
dE = 0.0
endif

end subroutine energy





