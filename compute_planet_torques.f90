subroutine compute_planet_torques
!
! For each planet, compute the torque it exerts on the disc
! cf Nayakshin (2014)
!

! Loop over each planet
do iplanet =1,nplanet


    do i = isr,ier

    ! Compute the Type II specific torque at this radius

    ! Compute the Type I specific torque

    ! Compute the interpolative factor f between migration regimes

    ! Compute the effective planet torque at this radius

    enddo


! Calculate the total torque exerted on the disc

total_planet_torque(:) = total_planet_torque(:) + torquei(:)

enddo


end subroutine compute_planet_torques