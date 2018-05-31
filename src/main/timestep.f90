subroutine timestep
!
!
! Determine maximum safe timestep for the simulation
! Must account for viscous evolution and planetary torques
! Timestep cannot exceed predefined maximum
!
!

use gravdata
use unitdata, only: yr,tdump

implicit none

integer :: i
real(kind=8) :: dtmin,C0,C1,dr2,dr


dtmin = 1.0d20
If (t .lt. 1.0d6) Then
   C0 = 0.001d0
Else If (t .lt. 1.0d7) Then
   C0 = 0.05d0
Else
   C0 = 0.25d0
EndIf
C1 = 0.5d0

do i = isr, ier
   dr  = (rf(i+1)-rf(i))
   dr2= dr**2
   dt  = C0 * dr2 / (6.0d0 * nu_tc(i))
   dtmin = min(dtmin,dt)
   dtmin = min(dtmin,0.01*tcool(i))
   dtmin = min(dtmin,maxstep*yr) ! timestep can not exceed maxstep
enddo

dt = dtmin
    
!print*, t/yr, dtmin/yr, 0.001*minval(pack(tcool, tcool>0.0))/yr
IF(dt>tdump) then
   print*, 'Reducing dt for rapid snapshotting'
   dt = tdump/2.0
endif

end subroutine timestep
