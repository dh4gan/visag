
!-----------------------------------------------------------
! Subroutine to write radial dumps 
!-----------------------------------------------------------

subroutine write_dump

  use gravdata
  use planetdata
  use unitdata
  use winddata

  implicit none

  real(kind=8) :: alpha_visc,sig_max,mdisk, tot_lumin,tcoolmin
  integer :: i,ifirst	

  write (*,103) 'Output at time t = ',t/yr
102 format (10E15.5)
103 format (A,E15.5)
104 format (4E15.5)

  ! Calculate parameters

  snew(:) = 0.0 ! need a dummy zeroed array to pass to calc_grav
  call calc_grav(snew)
  call luminosity(tot_lumin,spectrum,Tc,tau)

  ! write out disc profiles
  IF(counter<=9) THEN
     WRITE(fileno,'(I1)') counter
  ELSE IF(counter>9.and.counter<=99) THEN
     WRITE(fileno,'(I2)') counter
  ELSE IF(counter>99.and.counter<=999) THEN
     WRITE(fileno,'(I3)') counter
  ELSE IF(counter>999.and.counter<=9999) THEN
     WRITE(fileno,'(I4)') counter
  ELSE IF(counter>9999.and.counter<=99999) THEN
     WRITE(fileno,'(I5)') counter
  ENDIF
  fileno = TRIM(fileno)

  OPEN(iprof, file=TRIM(prefix)//'profile.'//fileno, status='unknown')
  write(iprof,*) nrgrid
  do i = isr, ier
     write (iprof,102) rz(i)/AU,sigma(i), cs(i), kappa(i), gamma(i),mu(i),Tc(i), tau(i), alpha_g(i),Q(i)
  enddo

  close(iprof)

  OPEN(ispec, file=TRIM(prefix)//'spectrum.'//fileno, status='unknown')
  write(ispec,*) 1000
  do i = 1, 1000
     write(ispec,*) spectrum(1,i), spectrum(2,i)
  enddo

  close(ispec)

  ! Compute disk mass and maximum surface density

  mdisk = 0.0d0
  sig_max = 0.0d0

  do i = isr, ier
     mdisk = mdisk + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma(i)
     sig_max = max(sig_max,sigma(i))
  enddo

  mdisk = mdisk / solarmass	
  tcoolmin = minval(tcool)
  ! Write out snapshot data (disk mass, sigma, total luminosity)
  write(itime,*) t/yr, dt, tcoolmin, mdisk, sig_max,tot_lumin

 
  return

end  subroutine write_dump
