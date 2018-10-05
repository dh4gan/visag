program semi_analytic_disc
  !-----------------------------
  ! Main program
  !
  ! Drives setup and evolution of viscous alpha disc
  ! 
  !----------------------------

  use gravdata
  use unitdata
  

  real(kind=8) :: tnext

  integer :: i, ifirst, ifirst2


  !	Setup grid, etc			
  call setup

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
        tnext = tdump*snapshotcounter			   	
     endif

  enddo

  close(itime)

end program semi_analytic_disc






