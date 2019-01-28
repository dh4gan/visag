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
  real :: dtorque, dTcdr, vr, dt_torque
  logical :: timestepOK

  integer :: i,j,k, iplanet
  integer :: Lcalc

  timestepOK = .false.

  ! Repeat calculation until timestep OK
  do while(.not.timestepOK)

     ! Determine disc parameters and viscosity

     snew(:) = 0.0
     Tnew(:) = 0.0

     call disc_properties

     if(planetchoice=='y') then
        ! Compute torques induced on the disc by planets
        call compute_planet_torques
     endif

     ! Get system timestep
     call timestep

     ! Set v_r = 0 outer bc
     ! This requires nu.sigma.r^{1/2} to have zero gradient

     nu_tc(isr-1) = nu_tc(isr)
     nu_tc(ier+1) = nu_tc(ier)
     sigma(ier+1) = sigma(ier) * sqrt(rz(ier)/rz(ier+1))

     ! Calculate net heating and cooling

     heatfunc(:) = 9.0*nu_tc(:)*sigma(:)*omegaK(:)*omegaK(:)/8.0


     ! Calculate wind parameters
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

           ! Symmetrised planet torque
           ! 1/2[ ((i+1) - (i)) + (i)-(i-1)) ]
           dtorque = 0.5*(torque_term(i+1) - torque_term(i-1))

           if(3.0*(term1-term2)-dtorque>0.0) then
              dtorque = 0.0
          endif
           
           do iplanet=1,nplanet
              if(i==iplanetrad(iplanet) .or. i-1==iplanetrad(iplanet))then
                 dtorque = 0.0
              endif
           enddo

        endif                

        vr = -3.0*(term2)/(rf(i)*sigma(i))

        dTcdr = (Tc(i+1) -Tc(i))*drzm1(i)   

        !if(t> 2.0*yr) then
        !write(75,'(9(1pe15.5,2X))') rz(i)/AU,sigma(i),total_planet_torque(i), &
        !     (total_planet_torque(i+1) - total_planet_torque(i))*drzm1(i), &
        !     torque_term(i), rzm1(i)*drzm1(i)*dtorque, &
        !     rzm1(i)*drzm1(i)*3.0*(term1-term2), &
        !     rzm1(i)*drzm1(i)*(3.0*(term1-term2) - dtorque)*dt, &
        !     nu_tc(i)             
        !endif

        snew(i) = sigma(i) + rzm1(i)*drzm1(i)*(3.0*(term1-term2) - dtorque)*dt -(sigdot_wind(i)-sigdot_accrete(i))*dt
        !write(*,'(3(1pe15.5,2X))') rz(i)/AU, 3.0*(term1-term2)-dtorque
        
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
    call move_planets
 endif

 
!print*, t/yr,dt/yr
!if(t>2.0*yr) STOP
!stop
  return
end subroutine evolve
