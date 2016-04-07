!-------------------------------------------------
! Subroutine to evolve the disk 
!-------------------------------------------------

subroutine evolve

  use gravdata
  use planetdata
  use unitdata


  real(kind=8) :: dtmin,C0,C1,dr2,dr
  real(kind=8) :: tot_lumin,term1,term2

  integer :: i,j,k
  integer :: Lcalc

  ! Determine disc parameters and viscosity

  snew(:) = 0.0
  Tnew(:) = 0.0

  call disc_properties

  ! Determine maximum safe timestep

  dtmin = 1.0d20
  If (t .lt. 1.0d6) Then
     C0 = 0.005d0
  Else If (t .lt. 1.0d7) Then
     C0 = 0.05d0
  Else
     C0 = 0.25d0
  EndIf
  C1    = 0.5d0

  do i = isr, ier
     dr  = (rf(i+1)-rf(i))
     dr2= dr**2
     dt  = C0 * dr2 / (6.0d0 * nu_tc(i))
     dtmin = min(dtmin,dt)
     dtmin = min(dtmin,0.01*tcool(i))
  enddo

  dt = dtmin
    
  IF(dt>tdump) dt = tdump/2.0

  ! Set v_r = 0 outer bc
  ! This requires nu.sigma.r^{1/2} to have zero gradient

  nu_tc(isr-1) = nu_tc(isr)
  nu_tc(ier+1) = nu_tc(ier)
  sigma(ier+1) = sigma(ier) * sqrt(rz(ier)/rz(ier+1))

  ! Calculate net heating and cooling

  heatfunc(:) = 9.0*nu_tc(:)*sigma(:)*omegaK(:)*omegaK(:)/8.0

  ! Evolve the surface density with the diffusive term

  !$OMP PARALLEL &
  !$OMP shared(rf1_2,nu_tc,sigma,rz1_2,dt,sigdot) &
  !$OMP shared(drzm1,drfm1,isr,ier) &
  !$OMP private(term1,term2,i)
  !$OMP DO SCHEDULE(runtime)
  do i = isr, ier
     
     term1=rf1_2(i+1)*(nu_tc(i+1)*sigma(i+1)*rz1_2(i+1)-nu_tc(i)*sigma(i)*rz1_2(i))*drfm1(i+1)
     term2=rf1_2(i)*(nu_tc(i)*sigma(i)*rz1_2(i)-nu_tc(i-1)*sigma(i-1)*rz1_2(i-1))*drfm1(i)

     vr = -3.0*term2/(rf(i)*sigma(i))
     dTcdr = (Tc(i+1) -Tc(i))*drzm1(i)

     snew(i) = sigma(i) + 3.0d0*rzm1(i)*drzm1(i)*(term1-term2)*dt +sigdot(i)*dt
     Tnew(i) = Tc(i) + 2.0*dt*(heatfunc(i)-coolfunc(i))/(cp(i)*sigma(i)) -vr*dTcdr*dt
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Copy back to surface density

  !$OMP PARALLEL &
  !$OMP shared(sigma,snew,isr,ier) &
  !$OMP private(i)
  !$OMP DO SCHEDULE(runtime)
  do i = isr, ier
     sigma(i) = snew(i)  
     Tc(i) = Tnew(i)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine evolve
