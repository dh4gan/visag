subroutine nbody_output(t)
! Outputs data to file
! Currently writes each particle to separate file

  use nbodydata
  use planetdata
use unitdata,only: yr,twopi

implicit none
real, intent(in) :: t
integer :: ibody,iplanet

102 format (1P,24E15.5)
103 format (1P, 7E15.5)

call nbody_acceleration(pos,vel,acc)
call nbody_system_properties

! output individual bodies to separate files

   do ibody=2,nbodies
     iplanet = ibody-1

     if(alive(iplanet)==0) cycle ! Skip finished particles

      write(ibody+inbodylog,102) t, mass(ibody),pos(:,ibody), vel(:,ibody), &
           acc(:,ibody),a(ibody), ecc(ibody), &
           inc(ibody), longascend(ibody), &
           argper(ibody), trueanom(ibody), &
           ekin(ibody),epot(ibody), etot(ibody),angmom(:,ibody), &
           tmig(iplanet)
      call flush(ibody)
   enddo

! Write log file containing global simulation data

write(inbodylog,103) t, dt_nbody/twopi, 100.0*maxerror/tolerance, system_energy, dE, system_ang, dL


end subroutine nbody_output
