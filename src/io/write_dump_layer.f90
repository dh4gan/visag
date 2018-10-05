!-----------------------------------------------------------
! Subroutine to write snapshots of the disc  
! as a function of radius
!
!-----------------------------------------------------------

subroutine write_dump

  use gravdata
  use magdata
  use unitdata

  implicit none

  real ::  alpha_visc,sig_max,mdisk, tot_lumin
  real :: mdot_grav, mdot_mag, grav_max,mag_max, mmag,mgrav 
  integer :: i,ifirst	

  write (*,103) 'Output at time t = ',t/yr
102 format (9E15.5)
103 format (A,E15.5)
104 format (4E15.5)
107 format(7E15.5)
110 format(10E15.5)
111 format (11E15.5)

  ! Calculate parameters

  call calc_layers
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

     write (iprof,110) rz(i)/AU,sigma(i), cs(i), kappa(i), gamma(i),mu(i),Tc(i), &
	tau(i),nu_tc(i), alpha_g(i)
  enddo

  close(iprof)

  OPEN(iprof,file=TRIM(prefix)//'layer.'//fileno,status='unknown')
  WRITE(iprof,*) nrgrid
  do i=isr,ier
     write(iprof,111) rz(i)/AU, sigma(i), sigma_m(i), sigma_tot(i), cs_m(i),kappa_m(i),gamma_m(i), &
	mu_m(i),tau_m(i),nu_m(i),alpha_m
  enddo
close(iprof)

  OPEN(ispec, file=TRIM(prefix)//'spectrum.'//fileno, status='unknown')
  write(ispec,*) 1000
  do i = 1, 1000
     write(ispec,*) spectrum(1,i), spectrum(2,i)
  enddo

  close(ispec)

  ! Compute disk mass and maximum surface density
  ! Also compute radially averaged accretion rates

  mdisk = 0.0d0
  mgrav= 0.0d0
  mmag = 0.0d0

  grav_max = 0.0d0
  mag_max = 0.0d0
  sig_max = 0.0d0

  do i = isr, ier
     mmag = mmag + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma_m(i)
     mgrav = mgrav + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma(i)
     mdisk = mdisk+2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma_tot(i)

     grav_max = max(grav_max,sigma(i))
     mag_max = max(mag_max,sigma_m(i))
     sig_max = max(sig_max,sigma_tot(i))

     mdot_grav = mdot_grav + 3.0*pi*nu_tc(i)*sigma(i)/REAL(nrgrid)
     mdot_mag = mdot_mag + 3.0*pi*nu_m(i)*sigma_m(i)/REAL(nrgrid)

  enddo

  mdisk = mdisk/solarmass	
  mgrav = mgrav/solarmass
  mmag = mmag/solarmass

  mdot_grav = mdot_grav*yr/solarmass
  mdot_mag = mdot_mag*yr/solarmass

  ! Write out snapshot data (disk mass, sigma, total luminosity)
  write(itime,110) t/yr, mgrav,mmag,mdisk, grav_max,mag_max, sig_max,tot_lumin, mdot_grav,mdot_mag

   open(87,file=TRIM(prefix)//'mri.'//fileno,status='unknown')
  do i=isr,ier
  write(87,*) rz(i)/AU, fullmri(i)
  enddo
  close(87)

  return

end  subroutine write_dump
