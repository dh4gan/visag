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
  real :: tmig1, lambda_dash



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

        ! Must first compute the total torque and then normalise
        ! to obtain the correct Type I migration timescale
           
        ! Integrate the torque over all radii
        typeInorm = typeInorm + exp(-deltap/(H+rhill))*sigma(i)/drzm1(i)        

     enddo

     ! Compute Type I migration timescale
     call calc_typeI_migration(iplanet,tmig1)

     ! Find normalisation constant to ensure correct Type I timescale

     !lambda_dash = mratio*omegaK(iplanetrad(iplanet))*ap(iplanet)*ap(iplanet)/(4.0*pi*G*tmig1*typeInorm)

     lambda_dash = ap(iplanet)*ap(iplanet)*omegaK(iplanetrad(iplanet))/(4.0*pi*G*tmig1*typeInorm)

     ! Now compute functional form of lambda 
     ! (assuming concentration of torque in planet local vicinity)
     do i=isr,ier
        deltap = abs(rz(i)-ap(iplanet))
        H = cs(i)/omegaK(i)  
        if(H>deltap) deltap =H
        if(rhill>deltap) deltap=rhill

        lambdaI(iplanet,i) = lambda_dash*exp(-deltap/(H+rhill))        
       
     enddo


     !**************************************************
     ! Now compute the relative dominance of each torque 
     !**************************************************

     do i= isr,ier

        !*****************************************************
        ! Compute the interpolative factor f between migration regimes
        !******************************************************

        ! Pressure criterion for Type I/ Type II

        Pcrit = 0.75 * H/rhill + 50.0*alpha_g(i)*(H/ap(iplanet))**2/mratio

        fII(iplanet,i) = exp(-Pcrit+1.0)

        if(fII(iplanet,i) > 1.0) fII(iplanet,i)=1.0

        !fII(iplanet,i) = 1.0 ! DEBUG LINE - REMOVE!
        
        !print*, i,fII(iplanet,i), lambdaI(iplanet,i), lambdaII(iplanet,i)
        !********************************************************
        ! Compute the total effective planet torque at this radius
        !*********************************************************

        torquei(iplanet,i) = lambdaI(iplanet,i)*(1.0-fII(iplanet,i)) + lambdaII(iplanet,i)*fII(iplanet,i)

     enddo


     ! Add planet contribution to the total torque exerted on the disc

     total_planet_torque(:) = total_planet_torque(:) + torquei(iplanet,:)

  enddo
  
   torque_term(:) = 2.0*omegaK(:)*rz(:)*rz(:)*sigma(:)*total_planet_torque(:)

end subroutine compute_planet_torques
