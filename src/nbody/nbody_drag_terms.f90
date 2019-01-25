subroutine nbody_drag_terms(position,velocity,acceleration)
! Calculates approximation to 3 drag terms given a migration timescale
! Migration drag at timescale tmig
! Eccentricity damping on timescale tmig/10
! Inclincation damping on timescale tmig/10
! (cf Alibert et al 2013)

! Forces not calculated on body 1 (the star)
! Forces not calculated if embryo is finished

  use nbodydata
  use planetdata, only: tmig, alive
implicit none

integer :: ibody,ix
real,dimension(3,nbodies),intent(in) :: position,velocity
real,dimension(3,nbodies),intent(inout) :: acceleration
real,dimension(nbodies) :: vdotr,rmag

! Migration drag first

do ix=1,3
    do ibody=2,nbodies
    if(tmig(ibody-1)>small .and. alive(ibody-1)==1) then
       acceleration(ix,ibody) = acceleration(ix,ibody)- &
            velocity(ix,ibody)/(2.0*tmig(ibody-1))
    endif
    enddo
enddo

! Eccentricity damping
vdotr(:) =  0.0
do ix=1,3
   vdotr(:) = vdotr(:)+ position(ix,:)*velocity(ix,:)
enddo

rmag(:) = 0.0

rmag(:) = position(1,:)*position(1,:) + &
     position(2,:)*position(2,:) + &
     position(3,:)*position(3,:)

do ix=1,3

do ibody=2,nbodies

 if (tmig(ibody-1)*rmag(ibody)>small .and.  alive(ibody-1)==1) then
acceleration(ix,ibody) = acceleration(ix,ibody) - &
     2.0*vdotr(ibody)*position(ix,ibody)/(rmag(ibody)*rmag(ibody)*dampfac*tmig(ibody-1))
endif
enddo
enddo

! Inclination damping

do ibody=2,nbodies
!if(embryo(ibody-1)%finished==1) cycle
if(tmig(ibody-1)>small .and. alive(ibody-1)==1) then
    acceleration(3,ibody) = acceleration(3,ibody) - &
         2.0*vel(3,ibody)/(dampfac*tmig(ibody-1))
endif
end do

end subroutine nbody_drag_terms
