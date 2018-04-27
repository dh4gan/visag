
!-----------------------------------------------------------
! Subroutine to write radial dumps 
!-----------------------------------------------------------

subroutine write_dump

  use gravdata
  use magdata
  use planetdata
  use unitdata
  use winddata

  implicit none

  real(kind=8) :: sig_max,mdisk, tot_lumin,tcoolmin
  real :: mdot_grav, mdot_mag, grav_max,mag_max, mmag,mgrav
  integer :: i,ifirst, iplanet

  write (*,103) 'Output at time t = ',t/yr
102 format (11(1PE15.5))
103 format (A,1PE15.5)
104 format (4(1PE15.5))
110 format(10(1PE15.5))
111 format (11(1PE15.5))


  ! Calculate disc properties and spectrum

  call disc_properties
  call luminosity(tot_lumin,spectrum,Tc,tau)

  write(fileno, snapshotformat)snapshotcounter

  fileno = TRIM(fileno)
! write out disc profiles
  OPEN(iprof, file=TRIM(prefix)//'_profile.'//fileno, status='unknown')
  write(iprof,*) t/yr, nrgrid
  do i = isr, ier
     write (iprof,102) rz(i)/AU,sigma(i), cs(i), kappa(i), gamma(i),mu(i),Tc(i), &
    tau(i), nu_tc(i), alpha_g(i),Q(i)
  enddo

  close(iprof)

  OPEN(ispec, file=TRIM(prefix)//'_spectrum.'//fileno, status='unknown')
  write(ispec,*) 1000
  do i = 1, 1000
     write(ispec,*) spectrum(1,i), spectrum(2,i)
  enddo

  close(ispec)

if (layerchoice=='y') then
OPEN(iprof,file=TRIM(prefix)//'layer.'//fileno,status='unknown')
WRITE(iprof,*) t/yr,nrgrid
do i=isr,ier
write(iprof,111) rz(i)/AU, sigma(i), sigma_m(i), sigma_tot(i), cs_m(i),&
    kappa_m(i),gamma_m(i), mu_m(i),tau_m(i),nu_m(i),alpha_m
enddo
close(iprof)
endif


if(planetchoice=='y') then
open(iprof, file=TRIM(prefix)//'_planets.'//fileno,status='unknown')
write(iprof,*)t/yr,nplanet,nactive
do iplanet=1,nplanet

    write(iprof,*) alive(iplanet),mp(iplanet)/mjup, &
         ap(iplanet)/AU, tmig(iplanet)/yr
enddo
close(iprof)
endif


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

mdot_wind = mdot_wind + sigdot_wind(i)*twopi*rz(i)/drzm1(i)

enddo

mdisk = mdisk/solarmass
mgrav = mgrav/solarmass
mmag = mmag/solarmass

mdot_grav = mdot_grav/msolyr
mdot_mag = mdot_mag/msolyr
mdot_wind = mdot_wind/msolyr

! Write out snapshot data (disk mass, sigma, total luminosity)
write(itime,111) t/yr, dt/yr, mdisk, tot_lumin, sig_max,mgrav,mmag,grav_max,mag_max,mdot_grav,mdot_mag,mdot_wind
call flush(itime)
 
  return

end  subroutine write_dump
