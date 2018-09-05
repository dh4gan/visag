SUBROUTINE layer_properties
  ! Subroutine calculates evolution of an MRI activated layer
  ! Checks for MRI activation at all disc radii
  ! either fully MRI to the midplane, or only the upper layer

  use gravdata
  use magdata
  use unitdata

  implicit none

  real(kind=8) :: coolfunclayer, twoDint,fine,oldtry
  real(kind=8) :: Teff,temp1,temp2,sig1,sig2
  real(kind=8) :: rho,Tnewlayer,heat,cool

  real(kind=8) :: try,T_try,rho_m, alpha_max,tcool_target

  real, parameter :: tolerance = 1.0d-2

  integer :: i, j, ntries

  !$OMP PARALLEL & 
  !$OMP shared(isr,ier,mstar) &
  !$OMP shared(sigma_tot,sigma,sigma_m) &
  !$OMP shared(cs_m,mu_m,nu_m,kappa_m,gamma_m) &
  !$OMP shared(Tc,rz,tau_m,T_source,tcool,omegaK) &
  !$OMP private(i,rho_m,T_try,try) &
  !$OMP private(temp1,temp2,sig1,sig2)
  !$OMP DO SCHEDULE(runtime)
  do i=isr,ier

     ! Calculate sigma_m, column density of magnetic layer

     ! If T drops below recombination temperature, MRI is turned off for fully MRI zones
     IF(Tc(i) < Toff) THEN
        fullmri(i)=fullmri(i)-0.1
        IF(fullmri(i) < 0.0) fullmri(i) =0.0
     ENDIF
     
     fullmri(i) = 1.0
     ! If T> Tcrit or MRI still active from previous step, GI layer disappears

     IF(Tc(i)> Tcrit.or.fullmri(i)/=0.0) THEN
        sigma_m(i) = sigma_tot(i)
        sigma(i) = 0.0
        alpha_g(i) = 0.0
        rho = 0.0
        IF(fullmri(i)==0.0) fullmri(i)=1.0

        mu_m(i) = 2.4
        kappa_m(i) = kappa(i)
        gamma_m(i) = 1.666
        mu_m(i) = 2.4

        ! Calculate equilibrium temperature for this layer by iteration
        ! Match heating and cooling from pure MRI turbulence

        ntries = 0

        try = 100.0
        temp1 = Tcrit
        temp2 = 1000.0*Tcrit

        DO WHILE (abs(try) > tolerance)

           ! Calculate state variables at this temperature
           cs_m(i) = SQRT(gamma_m(i)*k_B*Tc(i)/(mu_m(i)*m_H))
           rho_m = 0.5*sigma_m(i)*omegaK(i)/cs_m(i)

           call eos(rho_m,cs_m(i),T_try,kappa_m(i),mu_m(i),gamma_m(i))

           tau_m(i) = sigma_m(i)*kappa_m(i)

           ! Calculate midplane temperature derived from these variables
           nu_m(i) = alpha_m*cs_m(i)*cs_m(i)/omegaK(i)

           heat = 9.0*nu_m(i)*omegaK(i)*omegaK(i)*sigma_m(i)/8.0        
           Tnewlayer = heat*(tau_m(i) + 1.0/tau_m(i))/stefan
           Teff = Tc(i)**4/(tau_m(i)+1.0/(tau_m(i)))

           cool = stefan*Teff
           Tnewlayer = Tnewlayer**0.25
           Teff = Teff**0.25d0

           IF(heat > cool) temp1 = Tc(i)
           IF(heat < cool) temp2 = Tc(i)
           
           try = (Tnewlayer-Tc(i))/Tc(i)

           Tc(i) = (temp2+temp1)/2.0

           ntries = ntries + 1
           IF(ntries > 100 .and. (temp2-temp1)/temp1 < 1.0e-5) exit

        ENDDO

     ELSE IF(fullmri(i)==0) THEN

        ! Otherwise calculate MRI layer explicitly

        try = 100.0
        kappa_m(i) = 0.0
        gamma_m(i) = 1.666
        mu_m(i) = 2.4
        sigma_m(i) = 1.0d-2*sigma_tot(i)

        sig1 = sigma_m(i)/5.0
        sig2 = 5.0*sigma_m(i)

        ntries = 0
        fine = 0.1
        DO WHILE(abs(try).gt.tolerance)

           cs_m(i) = SQRT(gamma_m(i)*k_B*T_source(i)/(mu_m(i)*m_H))
           rho_m = 0.5*sigma_m(i)*omegaK(i)/cs_m(i)

           call eos(rho_m,cs_m(i),T_try,kappa_m(i),mu_m(i),gamma_m(i))

           tau_m(i) = sigma_m(i)*kappa_m(i)

           IF(tau_m(i) < tau_crit) sig1 = sigma_m(i)
           IF(tau_m(i) > tau_crit) sig2 = sigma_m(i)

           try = (sig2-sig1)/(sigma_m(i))

           sigma_m(i) = (sig1+sig2)/2.0

           ntries = ntries+1                

           IF(sigma_m(i)<=0.0) exit
        ENDDO

        sigma(i) = sigma_tot(i) - sigma_m(i)
        IF(sigma(i)<=0.0) sigma(i) = 0.0

     ENDIF
     nu_m(i) = alpha_m*cs_m(i)*cs_m(i)/omegaK(i)

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

return

END SUBROUTINE layer_properties
