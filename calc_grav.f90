  SUBROUTINE calc_grav(tauplus)
! Routine calculates state properties of the self-gravitating disc (layer)
! tauplus is passed in case there is extra optical depth above the layer

  use gravdata
  use unitdata

    implicit none
    
    real,dimension(nmax) :: tauplus
    real(kind=8) :: alpha_visc,  twoDint,fine,oldtry
    real(kind=8) :: rho,H, Teff
    integer :: i

    coolfunc(:) = 0.0
!$OMP PARALLEL &
  !$OMP shared(isr,ier,cs,sigma,kappa,mu,gamma,omegaK)&
  !$OMP shared(Tc,tau,nu_tc,tcool,T_source,alpha_g)&
  !$OMP shared(tauplus) &
  !$OMP private(H,rho,coolfunc,Teff,twoDint,i)
  !$OMP DO SCHEDULE(runtime)
  do i = isr, ier

     omegaK(i) = Sqrt(G*mstar/rz(i)**3.0d0)

     IF(Tc(i) < T_source(i)) Tc(i) = T_source(i)
     cs(i) = k_B*Tc(i)/(mu(i)*m_H)

     cs(i) = sqrt(cs(i))
     H = cs(i)/omegaK(i)

     IF(H/=0.0) THEN
        rho = 0.5d0*sigma(i)/H
     ELSE
        rho = 0.0
     ENDIF


     !	Interpolate over rho,T to get cs,kappa, mu, gamma

     IF(rho>=1.0e-25) THEN

        call eos_T(rho,cs(i),Tc(i),kappa(i),mu(i),gamma(i),cp(i))
        
        ! Calculate Q
        Q(i) = cs(i)*omegaK(i)/(pi*G*sigma(i))
!        write(77,*) t,i, rho,cs(i),H/AU, Tc(i),Q(i)

        ! Now calculate cooling time
        tau(i) = H*kappa(i)*rho + tauplus(i)

        ! Cooling function for disc
        coolfunc(i) = 16.0d0/3.0d0*stefan*(Tc(i)**4.0d0-T_source(i)**4.0d0)
        coolfunc(i) = coolfunc(i)*tau(i)/(1.0d0+tau(i)**2.0d0)

        Teff = Tc(i)**4.0d0*tau(i)/(1.0d0+tau(i)**2.0d0)
        Teff = Teff**0.25d0

        twoDint = cs(i)**2.0*sigma(i)/gamma(i)/(gamma(i)-1.0d0)

        If (coolfunc(i).ne.0.0d0) Then
           tcool(i) = twoDint/coolfunc(i)
        Else
           tcool(i) = 1.0d35
        EndIf

     ELSE
        Tc(i) = 0.0
        cs(i) = 0.0
        kappa(i) = 0.0
        mu(i) = 0.0
        gamma(i) = 0.0
        tcool(i) = 1.0d35
        coolfunc(i) = 0.0
     ENDIF

     !tcool(i) = tcool(i)*omegaK(i)

     alpha_g(i) = 0.01


    !IF(gamma(i)/=0.0) THEN
    !    alpha_g(i) =  9.0d0/4.0d0*gamma(i)*(gamma(i)-1.0)*tcool(i)*omegaK(i)
    !    alpha_g(i) =  1.0d0/alpha_g(i)
     
        !IF(t/yr> 999.0) print*, i, ' alpha is ', alpha_visc
    ! else if (alpha_g(i) .lt. 1.0d-12.or.gamma(i)==0.0) THEN
    !    alpha_g(i) = 1.0d-12
        
    ! endif

    ! IF(Q(i) < Qcrit) THEN
    !    alpha_g(i) = alpha_g(i)*(Qcrit*Qcrit/(Q(i)*Q(i)) - 1.0)
    ! ELSE
    !    alpha_g(i) = 0.0
    ! ENDIF

     nu_tc(i) = alpha_g(i)*cs(i)*cs(i)/omegaK(i)         

!     IF(t/yr > 900.0) print*, 'nu_tc ',i,t/yr,dt/yr,nu_tc(i),alpha_visc,tcool,rho,cs(i),Tc(i),sigma(i),sigma_m(i)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

return

END SUBROUTINE calc_grav
