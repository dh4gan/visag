!-------------------------------------------------
! Subroutine to evolve a two layer disk
!-------------------------------------------------

subroutine evolve_layers

  use gravdata
  use magdata
  use unitdata


  implicit none

  real(kind=8) :: dtmin,C0,C1,dr2,dr,mdisk,sig_max
  real(kind=8) :: tot_lumin,term1,term2,nu_try
  integer :: i,j,k, imin
  integer :: Lcalc

  ! Determine disc properties

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
     !dt  = C0 * dr2 / (6.0d0 * nu_tc(i))
     !dtmin = min(dtmin,dt)
     nu_try = max(nu_tc(i),nu_m(i))
     IF(nu_try/=0.0) THEN
        dt = C0*dr2/(6.0d0*nu_try)
        dtmin = min(dtmin,dt)
        imin = i
     ENDIF
  enddo

  dt = dtmin

  IF(dt>tdump) dt = tdump/2.0

  ! Set v_r = 0 outer bc
  ! This requires (nu_g sigma_g + nu_m sigma_m) r^{1/2} to have zero gradient
  ! 
  nu_tc(isr-1) = nu_tc(isr)
  nu_tc(ier+1) = nu_tc(ier)

  nu_m(isr-1) = nu_m(isr)
  nu_m(ier+1) = nu_m(ier)

  ! Require constraints on sigma_m to close system
  ! Assume v_r due to sigma only is also zero

  sigma(ier+1) = sigma(ier)* sqrt(rz(ier)/rz(ier+1))
  sigma_m(ier+1) = 0.0

  sigma_tot(ier+1) = sigma(ier+1) + sigma_m(ier+1)

  ! Evolve the surface density with the diffusive term


  !$OMP PARALLEL &
  !$OMP shared(rf1_2,nu_tc,sigma,rz1_2,dt,sigdot) &
  !$OMP shared(sigma_m, nu_m) &
  !$OMP shared(drzm1,drfm1,isr,ier) &
  !$OMP private(term1,term2,i)
  !$OMP DO SCHEDULE(runtime)
  do i = isr, ier
     term1=rf1_2(i+1)*((nu_tc(i+1)*sigma(i+1)+nu_m(i+1)*sigma_m(i+1))*rz1_2(i+1)- &
          (nu_tc(i)*sigma(i)+nu_m(i)*sigma_m(i))*rz1_2(i))*drfm1(i+1)
     term2=rf1_2(i)*((nu_tc(i)*sigma(i)+nu_m(i)*sigma_m(i))*rz1_2(i)- & 
          (nu_tc(i-1)*sigma(i-1)+nu_m(i-1)*sigma_m(i-1))*rz1_2(i-1))*drfm1(i)
     sigma_tot(i) = sigma(i) +sigma_m(i)+ 3.0d0*rzm1(i)*drzm1(i)*(term1-term2)*dt+sigdot(i)*dt
        !IF(mag_switch==1.0.and.i==isr) print*, sigma_tot(i), sigma_m(i),sigma(i), term1-term2, &
        !     3.0d0*rzm1(i)*drzm1(i)*(term1-term2)*dt
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine evolve_layers
