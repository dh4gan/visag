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
  real :: rhill,H,mratio,deltap, Pcrit, typeInorm
  real :: deltap_next, deltap_prev
  real :: Hprev, Hnext, dr_sigmaH,tmig1


  torquei(:,:) = 0.0
  total_planet_torque(:) = 0.0

  ! Loop over each planet
  do iplanet =1,nplanet

     mratio = mp(iplanet)/Mstar
     rhill = ap(iplanet)*(mratio/3.0)**0.333


     typeInorm = 0.0    

     !**************************************************************
     ! Compute Type II and Type I torques (appropriately normalised)
     !**************************************************************

     do i = isr,ier


        !*****************************************************
        ! Compute the Type II specific torque at this radius
        !*****************************************************

        deltap = abs(rz(i)-ap(iplanet))
        H = cs(i)/omegaK(i)        

        if(H>deltap) deltap =H

        if(rhill>deltap) deltap=rhill

        lambdaII(iplanet,i) = 0.5*mratio*mratio/deltap**4
        if(rz(i) < ap(iplanet)) then
           lambdaII(iplanet,i) = lambdaII(iplanet,i)*rz(i)**4
        else
           lambdaII(iplanet,i) = lambdaII(iplanet,i)*ap(iplanet)**4
        endif

        !*************************************
        ! Compute the Type I specific torque (cf Nayakshin)
        !*************************************

        ! Type I torque depends on gradient of sigma*H^2

        if(i>isr .and. i <ier) then

           Hprev = cs(i-1)/omegaK(i-1)
           Hnext = cs(i+1)/omegaK(i+1)


           deltap_prev = abs(rz(i-1)-ap(iplanet))
           deltap_next = abs(rz(i+1)-ap(iplanet))
           Hprev = cs(i-1)/omegaK(i-1)
           Hnext = cs(i+1)/omegaK(i+1)

           dr_sigmaH = 0.5*drzm1(i)*(sigma(i+1)*Hnext*Hnext*exp(-deltap_next/(H+rhill)) - &
                2.0*sigma(i)*H*H*exp(-deltap/(H+rhill)) - &
                sigma(i-1)*Hprev*Hprev*exp(-deltap_prev/(H+rhill)))

        else
           dr_sigmaH = 0.0
        endif


        lambdaI(iplanet,i) = 0.25*mratio*dr_sigmaH/(sigma(i)*ap(iplanet)*ap(iplanet))

        typeInorm = typeInorm + exp(-deltap/(H+rhill))/drzm1(i)

        lambdaI(iplanet,i) = lambdaI(iplanet,i)*exp(-deltap/(H+rhill))

     enddo

     ! Normalise Type I torque

     lambdaI(iplanet,:) = lambdaI(iplanet,:)/typeInorm


     !**************************************************
     ! Now compute the relative dominance of each torque 
     !**************************************************

     do i= isr,ier

        !*****************************************************
        ! Compute the interpolative factor f between migration regimes
        !******************************************************

        ! Pressure criterion for Type I/ Type II

        Pcrit = 0.75 * H/rhill + 50.0*alpha_g(i)*(H/ap(iplanet))**2/mratio

        fII(iplanet,i) = exp(-Pcrit-1.0)
        if(fII(iplanet,i) > 1.0) fII(iplanet,i)=1.0


        !fII(iplanet,i) = 1.0 ! DEBUG LINE - REMOVE!

        !********************************************************
        ! Compute the total effective planet torque at this radius
        !*********************************************************

        torquei(iplanet,i) = lambdaI(iplanet,i)*(1.0-fII(iplanet,i)) + lambdaII(iplanet,i)*fII(iplanet,i)

     enddo


     ! Add planet contribution to the total torque exerted on the disc

     total_planet_torque(:) = total_planet_torque(:) + torquei(iplanet,:)

     torque_term(:) = 2.0*omegaK(:)*omegaK(:)*rz(:)*rz(:)*sigma(:)*total_planet_torque(:)

  enddo

end subroutine compute_planet_torques
