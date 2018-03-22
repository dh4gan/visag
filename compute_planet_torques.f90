subroutine compute_planet_torques
!
! For each planet, compute the torque it exerts on the disc
! cf Nayakshin (2014)
!

use planetdata
use gravdata
use unitdata

implicit none

integer :: iplanet,i
real :: rhill,H,mratio,deltap, Pcrit

torquei(:,:) = 0.0
total_planet_torque(:) = 0.0

! Loop over each planet
do iplanet =1,nplanet

    mratio = mp(iplanet)/Mstar
    rhill = ap(iplanet)*(mratio/3.0)**0.333

    do i = isr,ier

    !*****************************************************
    ! Compute the Type II specific torque at this radius
    !*****************************************************

    deltap = abs(rz(i)-ap(iplanet))
    H = cs(i)/omegaK(i)

    if(H>deltap) deltap =H

    lambdaII(iplanet,i) = 0.5*mratio*mratio/(deltap**4)

    if(rz(i) < ap(iplanet)) then
        lambdaII(iplanet,i) = lambdaII(iplanet,i)*rz(i)**4
    else
        lambdaII(iplanet,i) = lambdaII(iplanet,i)*ap(iplanet)**4
    endif

    !*************************************
    ! Compute the Type I specific torque
    !*************************************

    tmigI(iplanet) = 1.0 ! TODO compute type I migration timescale
    lambdaI(iplanet,i) = 0.5*sqrt(G*Mstar*ap(iplanet))/tmigI(iplanet)

    lambdaI(iplanet,i) = lambdaI(iplanet,i)*exp(-deltap/(H+rhill)) ! TODO is this deltap or deltar? What is deltar?

    !*****************************************************
    ! Compute the interpolative factor f between migration regimes
    !******************************************************

    ! Pressure criterion for Type I/ Type II

    Pcrit = 0.75 * H/rhill + 50.0*alpha_g(i)*(H/ap(iplanet))**2/mratio

    fII(iplanet,i) = exp(1.0-Pcrit)
    if(fII(iplanet,i) < 1.0) fII(iplanet,i)=1.0

    !********************************************************
    ! Compute the total effective planet torque at this radius
    !*********************************************************

    torquei(iplanet,i) = lambdaI(iplanet,i)*(1.0-fII(iplanet,i)) + lambdaII(iplanet,i)*fII(iplanet,i)

    enddo


! Add planet contribution to the total torque exerted on the disc

total_planet_torque(:) = total_planet_torque(:) + torquei(iplanet,:)

enddo


end subroutine compute_planet_torques