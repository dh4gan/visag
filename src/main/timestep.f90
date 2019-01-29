subroutine timestep
!
!
! Determine maximum safe timestep for the simulation
! Must account for viscous evolution and planetary torques
! Timestep cannot exceed predefined maximum
!
!

use gravdata
use unitdata, only: yr,twopi, tdump, nbodychoice
use planetdata, only: torque_term, planetchoice
use nbodydata, only: dt_nbody, pos,vel

implicit none

integer :: i
real(kind=8) :: dtmin,C0,C1,dr2,dr
real :: dtvisc, dtmin_visc, dttorq, dtmin_torque


dtmin_visc = 1.0e30
dtmin_torque = dtmin_visc

If (t .lt. 1.0d6) Then
   C0 = 0.01d0
Else If (t .lt. 1.0d7) Then
   C0 = 0.05d0
Else
   C0 = 0.25d0
EndIf
C1 = 0.5d0

!C0 = C0/10.0

do i = isr, ier

   ! Compute minimum viscous timestep
   dr  = (rf(i+1)-rf(i))
   dr2 = dr**2
   dtvisc  = C0 * dr2 / (6.0d0 * nu_tc(i))
   dtmin_visc = min(dtmin_visc,dtvisc)
   !dtmin = min(dtmin,0.001*tcool(i))

   ! If planets present, compute minimum timestep due to torques
   if(planetchoice=='y') then
      if(torque_term(i)>1.0e-30) then
         dttorq = C0*dr*rz(i)*sigma(i)/torque_term(i)
      else
         dttorq = 1.0e30
      endif

      dtmin_torque = min(dtmin_torque, dttorq)   
   endif

enddo

dt = dtmin_visc
dt = min(dtmin_visc, dtmin_torque)

! If doing an N Body run, check against minimum N Body timestep
if(nbodychoice=="y") then

   call nbody_timestep(pos,vel)
   dt = min(dt, dt_nbody*yr/twopi)  
endif

dt = min(dt,maxstep*yr) ! timestep can not exceed maxstep

if(dt>tdump) then
   print*, 'Reducing dt for rapid snapshotting'
   dt = tdump/10.0
endif

print*, dt/yr, maxstep, dtmin_torque/yr, dtmin_visc/yr, dt_nbody/twopi

! Make sure N Body timestep is corrected
if(nbodychoice=="y") dt_nbody = dt*twopi/yr

end subroutine timestep
