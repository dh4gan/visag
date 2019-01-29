program semi_analytic_disc
  !-----------------------------
  ! Main program
  !
  ! Drives setup and evolution of viscous alpha disc
  ! 
  !----------------------------

  use gravdata
  use unitdata
  use nbodydata, only: pos,vel
  

  real(kind=8) :: tnext

  integer :: i, ifirst, ifirst2

  !	Setup grid, etc			
  call setup

  ! If running in N Body mode, initialise N Body timestep
  call nbody_timestep(pos,vel)

  !	Write first dump
  t = 0.0d0
  snapshotcounter = 1

  call write_dump

  print '(a)', 'Initial Conditions Written to File'

  tnext = tdump

  ! Begin the simulation
  do while (t .lt. trun)

    if(runmode=='l') then

     call evolve_layers ! IN DEVELOPMENT
  else
     call evolve
  endif

     t = t + dt           
  
     If (t .gt. tnext) then
        snapshotcounter = snapshotcounter +1
        call write_dump

        if(nbodychoice=="y") call nbody_output(t)

        
        tnext = tdump*snapshotcounter


        
     endif

  enddo

  call nbody_deallocate_arrays

  close(itime)

end program semi_analytic_disc






