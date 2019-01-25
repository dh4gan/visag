subroutine nbody_output(t)
! Outputs data to file
! Currently writes each particle to separate file

use embryodata
use eosdata,only: yr,twopi

implicit none
real, intent(in) :: t
integer :: ibody,iembryo

102 format (1P,24E15.5)
103 format (1P, 7E15.5)

call nbody_acceleration(pos,vel,acc)
call nbody_system_properties

! output individual bodies to separate files

   do ibody=2,nbodies
     iembryo = ibody-1

     if(embryo(iembryo)%finished==1) cycle ! Skip finished particles

      write(ibody+inbodylog,102) t/yr, mass(ibody),pos(:,ibody), vel(:,ibody), &
           acc(:,ibody),embryo(iembryo)%semimaj, embryo(iembryo)%ecc, &
           embryo(iembryo)%inc, embryo(iembryo)%longascend, &
           embryo(iembryo)%argper, embryo(iembryo)%trueanom, &
           ekin(ibody),epot(ibody), etot(ibody),angmom(:,ibody), &
           embryo(iembryo)%tmig
      call flush(ibody)
   enddo

! Write log file containing global simulation data

write(inbodylog,103) t/yr, dt_nbody/twopi, 100.0*maxerror/tolerance, system_energy, dE, system_ang, dL


end subroutine nbody_output
