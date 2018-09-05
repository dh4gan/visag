program semi_analytic_disc

  ! Code to determine surface density of alpha disc.

  use gravdata
  use unitdata
  

  real(kind=8) :: trem_dump

  integer :: i, ifirst, ifirst2


  !	Setup grid, etc			
  call setup

  !	Write first dump
  t = 0.0d0
  snapshotcounter = 1
  call write_dump

  print '(a)', 'Initial Conditions Written to File'

  trem_dump = tdump

  ! Begin the simulation
  do while (t .lt. trun)

    if(runmode=='l') then

     call evolve_layers
else
    call evolve
endif

     t = t + dt
       
     trem_dump = trem_dump - dt			
  
     If (trem_dump .lt. 0.0d0) then
        snapshotcounter = snapshotcounter +1
        call write_dump			
        trem_dump = tdump			   	
     endif

  enddo

  close(itime)

end program semi_analytic_disc






