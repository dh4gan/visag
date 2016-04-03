program semi_analytic_disc

  ! Code to determine surface density of alpha disc.

  use gravdata
  use unitdata
  

  real(kind=8) :: trem_dump

  integer :: i, ifirst, ifirst2
  logical :: disk_exist

  !	Check that parameter file exists

  inquire(file='disc.par',exist=disk_exist)

  IF(disk_exist.eqv..false.) THEN
     print*, 'ERROR! input file disc.par does not exist!'	
     STOP
  ENDIF

  !	Setup grid, etc			
  call setup

  !	Write first dump
  t = 0.0d0
  counter = 1
  call write_dump	
  print*, 'Initial Conditions Written to File'
  trem_dump = tdump

  do while (t .lt. trun)			
     call evolve
     t = t + dt
  
     trem_dump = trem_dump - dt			
  
     If (trem_dump .lt. 0.0d0) then
        counter = counter +1
        call write_dump			
        trem_dump = tdump			   	
     endif

  enddo

  close(itime)

end program semi_analytic_disc






