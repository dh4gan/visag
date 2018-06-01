!-------------------------------------------------
! Subroutine to evolve the disk 
!-------------------------------------------------

subroutine evolve

  use gravdata
  use planetdata
  use unitdata
use winddata, only: sigdot_wind, sigdot_accrete

  implicit none

  real(kind=8) :: dtmin,C0,C1,dr2,dr
  real(kind=8) :: tot_lumin,term1,term2
  real :: dtorque, dTcdr, vr
  logical :: timestepOK

  integer :: i,j,k
  integer :: Lcalc

  timestepOK = .false.

  ! Repeat calculation until timestep OK
  do while(.not.timestepOK)

     ! Determine disc parameters and viscosity

     snew(:) = 0.0
     Tnew(:) = 0.0

     call disc_properties

     ! Get system timestep
     call timestep

     ! Set v_r = 0 outer bc
     ! This requires nu.sigma.r^{1/2} to have zero gradient

     nu_tc(isr-1) = nu_tc(isr)
     nu_tc(ier+1) = nu_tc(ier)
     sigma(ier+1) = sigma(ier) * sqrt(rz(ier)/rz(ier+1))

     ! Calculate net heating and cooling

     heatfunc(:) = 9.0*nu_tc(:)*sigma(:)*omegaK(:)*omegaK(:)/8.0

     if(planetchoice=='y') then
        ! Compute torques induced on the disc by planets
        call compute_planet_torques
     endif


     call compute_wind

     ! Evolve the surface density with the diffusive term
     
     !$OMP PARALLEL &
     !$OMP shared(rf1_2,nu_tc,sigma,rz1_2,dt,sigdot) &
     !$OMP shared(drzm1,drfm1,isr,ier) &
     !$OMP private(term1,term2,i)
     !$OMP DO SCHEDULE(runtime)
     do i = isr, ier
        
        term1=rf1_2(i+1)*(nu_tc(i+1)*sigma(i+1)*rz1_2(i+1)-nu_tc(i)*sigma(i)*rz1_2(i))*drfm1(i+1)
        term2=rf1_2(i)*(nu_tc(i)*sigma(i)*rz1_2(i)-nu_tc(i-1)*sigma(i-1)*rz1_2(i-1))*drfm1(i)
        
        dtorque = 0.0
        if(planetchoice=='y') then
           dtorque = torque_term(i+1)- 2.0*torque_term(i) - torque_term(i)
        endif
        
        !print*, dtorque/(term1-term2)
        vr = -3.0*term2/(rf(i)*sigma(i))
        dTcdr = (Tc(i+1) -Tc(i))*drzm1(i)
        
        snew(i) = sigma(i) + rzm1(i)*drzm1(i)*(3.0*(term1-term2) - dtorque)*dt -(sigdot_wind(i)-sigdot_accrete(i))*dt
        
        if(snew(i)<0.0) snew(i) = 0.0
        if(runmode/='Q') then
           Tnew(i) = Tc(i) + 2.0*dt*(heatfunc(i)-coolfunc(i))/(cp(i)*sigma(i)) -vr*dTcdr*dt
        else
           Tnew(i) = Tc(i)
        endif
     
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     
     ! Check if surface density negative anywhere - if so, repeat calculation with reduced timestep
     timestepOK = .true.
     do i=isr,ier
        !print*, i, snew(i)
        if(snew(i)<0.0)  timestepOK = .false.

     enddo


     if(.not.timestepOK) then
        print*, 'Reducing maximum allowed timestep ',maxstep, maxstep/10.0
        maxstep = maxstep/10.0

        if(maxstep<1.0e-10) then
           print*, 'MAXIMUM TIMESTEP TOO LOW: ending run'
           stop
        endif
     endif

  enddo
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


   if(planetchoice=='y') then
    ! Move planets
    call migrate_planets
endif

  return
end subroutine evolve
