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
  real :: rhill,mratio,deltap, Pcrit, typeInorm
  real :: tmig1, lambda_dash, deltamax,aspectratio, softenfactor
  logical :: soften

  call find_planets_in_disc

  torquei(:,:) = 0.0
  total_planet_torque(:) = 0.0
  torque_term(:) = 0.0

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

        deltamax = H(i)

        soften = .false.
        if(1.44*rhill>deltamax) deltamax = 1.44*rhill
        
        if(deltap<deltamax) then
           deltap =deltamax
           soften = .true.
        endif

        lambdaII(iplanet,i) = 0.5*mratio*mratio/(deltap)**4
        if(rz(i) < ap(iplanet)) then
           lambdaII(iplanet,i) = -lambdaII(iplanet,i)*(rz(i))**4
        else
           lambdaII(iplanet,i) = lambdaII(iplanet,i)*(ap(iplanet))**4
        endif

        if(soften) then

           softenfactor = abs(rz(i)-ap(iplanet))/deltamax
           if(softenfactor < 0.0001) softenfactor = 0.0001           
           lambdaII(iplanet,i)=lambdaII(iplanet,i)*softenfactor           
        endif

        !*************************************
        ! Compute the Type I specific torque (cf Nayakshin)
        !*************************************

        ! Must first compute the total torque and then normalise
        ! to obtain the correct Type I migration timescale
           
        ! Integrate the torque over all radii

        if(rz(i)<ap(iplanet)) then
           typeInorm = typeInorm + exp(-deltap/(H(i)+rhill))*sigma(i)/drzm1(i)    
        else
           typeInorm = typeInorm + exp(-deltap/(H(i)+rhill))*sigma(i)/drzm1(i)    
        endif
       
     enddo

     ! Compute Type I migration timescale
     call calc_typeI_migration(iplanet,tmig1)

     ! Find normalisation constant to ensure correct Type I timescale

     if(tmig1 > 1.0e-30) then
        lambda_dash = ap(iplanet)*ap(iplanet)*omegaK(iplanetrad(iplanet))/(4.0*pi*G*tmig1*typeInorm)
      else
        lambda_dash = 0.0
      endif

     ! Now compute functional form of lambda 
     ! (assuming concentration of torque in planet local vicinity)
     do i=isr,ier

        deltap = abs(rz(i)-ap(iplanet))

        if(deltap<H(i)) deltap =H(i)
        if(deltap<rhill) deltap=rhill

        lambdaI(iplanet,i) = lambda_dash*exp(-deltap/(H(i)+rhill))

        !if(rz(i)<ap(iplanet)) lambdaI(iplanet,i) = -lambdaI(iplanet,i)
       
     enddo


     !**************************************************
     ! Now compute the relative dominance of each torque 
     !**************************************************

     Pcrit = 0.75 * H(iplanetrad(iplanet))/rhill + 50.0*alpha_g(iplanetrad(iplanet))*(H(iplanetrad(iplanet))/ap(iplanet))**2/mratio

     

     fII(iplanet) = exp(-Pcrit+1.0)
     if(fII(iplanet) > 1.0) fII(iplanet)=1.0

     fII(iplanet) = 1.0 ! DEBUG LINE - REMOVE!

      !print*, 'Desired migration rate: ', t/yr, tmig1/yr, Pcrit, fII(iplanet)
 
      !********************************************************
      ! Compute the total effective planet torque at this radius
      !*********************************************************

     torquei(iplanet,:) = lambdaI(iplanet,:)*(1.0-fII(iplanet)) + lambdaII(iplanet,:)*fII(iplanet)

     ! Add planet contribution to the total torque exerted on the disc

     total_planet_torque(:) = total_planet_torque(:) + torquei(iplanet,:)

  enddo

  ! Ensure torque magnitude does not exceed maximum permitted value
  ! (Alexander & Armitage 2009, ApJ, 704, 989)

  do i=isr,ier
     aspectratio = H(i)/rz(i)

     if(abs(total_planet_torque(i)) > 0.1*aspectratio) then
        total_planet_torque(i) = 0.1*aspectratio*total_planet_torque(i)/abs(total_planet_torque(i))

        !print*, rz(i)/AU, aspectratio, total_planet_torque(i)
     endif
  enddo
       

  ! Set brief time delay for planet torque activation

  if(t < tdelay_planettorque*yr) then
     total_planet_torque(:) = total_planet_torque(:)*1.0e-30
  endif

   torque_term(:) = 2.0*omegaK(:)*rf(:)*rf(:)*sigma(:)*total_planet_torque(:)

   ! Zero torque_term at planet locations

   do iplanet=1,nplanet
      torque_term(iplanetrad(iplanet)) = 0
      torque_term(iplanetrad(iplanet)+1) = 0
   enddo

end subroutine compute_planet_torques
