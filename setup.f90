
!------------------------------------------------
! Subroutine to setup initial conditions and grid
!------------------------------------------------

subroutine setup

  use gravdata
  use magdata
  use planetdata
  use winddata
  use unitdata

  implicit none

  real(kind=8) :: rmax,sig_r
  real(kind=8) :: mdisk,mdisktry,dz
  real(kind=8) :: beta,rwindout,area, sig_0_try,tolerance

  integer :: i,ngrid, nzeros

  logical :: disk_exist

  !	Check that parameter file exists

  inquire(file=paramfile,exist=disk_exist)

  if(disk_exist.eqv..false.) then
   print '(a,a,a)', 'ERROR! input file ',paramfile,' does not exist!'
   stop
  endif

  open (unit=10,file=paramfile,form='formatted',status='old')

  read(10,*) prefix  ! File output prefix
  read(10,'(a1)') runmode ! f = fixed alpha, g=self-gravitating, Q = self-gravitating, fixed Q
  read(10,'(a1)') layerchoice ! Run this with an MRI upper layer? (y/n)
  read(10,*) alpha_visc ! If fixed alpha viscosity, define it here
  read(10,*) trun    ! Maximum runtime
  read(10,*) tdump   ! Snapshot Interval
  read(10,*) nrgrid  ! Number of r cells
  read(10,*) nzgrid  ! Number of z cells
  read(10,*) mstar   ! Star mass
  read(10,*) mdisk   ! Disc Mass
  read(10,*) sig_r   ! Initial surface density profile
  read(10,*) T0      ! Temperature at 1 AU
  read(10,*) p_T     ! Temperature profile
  read(10,*) rin     ! Inner disc radius
  read(10,*) rout    ! Spectral Break
  read(10,*) rmax    ! Outer disc radius
  read(10,*) zmax    ! Maximum altitude
  read(10,*) mdot_init ! Initial accretion rate
  read(10,*) rwind    ! Wind Radius
  read(10,*) mdot_wind ! Wind loss rate
  read(10,*) rremove  ! Radius at which planets are removed

  prefix = trim(prefix)

  write(*,*) " "
  write(*,*) "-----------------------------------------------"
if(runmode=='g') then
  write(*,*) "           VISCOUS SELF-GRAVITATING DISC CODE"
else if(runmode=='Q')then
    write(*,*) "         FIXED-Q VISCOUS SELF-GRAVITATING DISC CODE"
else
    write(*,*) "         FIXED-ALPHA VISCOUS DISC CODE"
endif
  write(*,*) "     Modified by D.Forgan, 29th April 2010       "
  write(*,*) "-----------------------------------------------"
  write(*,*) " "
  write(*,*) "-----------------------------------------------"
  write(*,*) " Input file: ./",trim(paramfile)
  write(*,*) "-----------------------------------------------"
  write(*,*) " "
  write(*,*) " "			
  write(*,*) "-----------------------------------------------"

  if(layerchoice=='y') then
    write(*,*) 'Disc will run with an MRI active upper layer'
endif
  write (*,100) ' - nrgrid  = ',nrgrid
  write(*,100)  ' - nzgrid  = ',nzgrid
  write(*,101) ' - mstar = ', mstar, ' solar masses'
  write(*,101) ' initial mdisk = ', mdisk, ' solar masses'
  write(*,101) ' Disc sigma power law = ',sig_r
  write (*,101) ' - rin    = ',rin, ' AU'
  write (*,101) ' - rout   = ',rout, ' AU'
  write(*,101) ' rmax(profile break) = ',rmax, ' AU'
  write(*,101)  ' - zmax   = ',zmax, ' AU'
  write (*,101) ' - trun   = ',trun, ' yr'
  write (*,101) ' - tdump  = ',tdump, ' yr'	
  write (*,101) ' - mdot   = ',mdot_init, ' solar masses / yr'
  write (*,101) ' - r_wind = ',rwind, ' AU'
  write (*,101) ' - wind loss = ',mdot_wind, ' solar masses / yr'
  write (*,101) ' - inner cut = ',rremove, ' AU'
  write(*,*) "-----------------------------------------------"
100 format (A,I5)
101 format (A,E14.3,A)


  ! Predicted number of files, and resulting format

    nfiles = trun/tdump
    nzeros = int(log10(nfiles)) +2
    write(zerostring, '(I1)')nzeros
    snapshotformat = "(I"//TRIM(zerostring)//"."//TRIM(zerostring)//")"

  ! Convert variables to cgs units

  mdot_init = mdot_init * msolyr
  mdot_wind = mdot_wind * msolyr

  rin = rin*AU
  rout = rout*AU
  rwind = rwind*AU
  rmax = rmax *AU
  zmax = zmax*AU

  trun = trun*yr
  tdump = tdump*yr

  mstar = mstar*solarmass
  mdisk = mdisk*solarmass

! Allocate arrays

  nmax = nrgrid + 10

allocate(sigma(nmax))
allocate(snew(nmax))
allocate(Tc(nmax))
allocate(Tnew(nmax))
allocate(cs(nmax))
allocate(tau(nmax))
allocate(nu_tc(nmax))
allocate(gamma(nmax))
allocate(mu(nmax))
allocate(kappa(nmax))
allocate(tcool(nmax))
allocate(alpha_g(nmax))
allocate(omegaK(nmax))
allocate(T_source(nmax))
allocate(Q(nmax))
allocate(heatfunc(nmax))
allocate(coolfunc(nmax))
allocate(cp(nmax))
allocate(sigdot(nmax))

! EDIT - Watch this !
cp(:) = 0.0
mu(:) = 2.4
gamma(:) = 7.0/5.0

! Set up grid arrays

allocate(rz(nmax))
allocate(rzm1(nmax))
allocate(rz1_2(nmax))
allocate(drzm1(nmax))

allocate(rf(nmax))
allocate(rf1_2(nmax))
allocate(drfm1(nmax))
allocate(zgrid(nmax))


! Set up layer arrays

allocate(sigma_m(nmax))
allocate(sigma_tot(nmax))
allocate(cs_m(nmax))
allocate(tau_m(nmax))
allocate(nu_m(nmax))
allocate(gamma_m(nmax))
allocate(mu_m(nmax))
allocate(kappa_m(nmax))
allocate(fullmri(nmax))

sigma_m(:) = 0.0
sigma_tot(:) = 0.0
mag_switch = 0.0
fullmri(:) = 0.0

  sigma(:) = 0.0
  Tc(:) = 0.0
  sigdot(:) = 0.0

! Set up source T according to input data
! TODO - add other temperature options

  T_source(:) = 10.0

  !	Read in EOS
  call eosread

  open(itime,file=TRIM(prefix)//'.log',status='unknown')

  ! Set up the basic grid
  ! rz denotes zone centered radius
  ! rf denotes face centered, displaced 1/2 cell inwards
  ! also precompute some helpful quantitites

  beta = (rout/rin)**(1.0d0/dble(nrgrid))
  dz = zmax/dble(nzgrid)

  isr = 2
  ier = nrgrid+1

  isz = 1
  iez = nzgrid

  rf(isr) = rin
  rz(isr) = rin * sqrt(beta)
  rzm1(isr)  = 1.0d0 / rz(isr)
  rf1_2(isr) = sqrt(rf(isr))
  rz1_2(isr) = sqrt(rz(isr))

  rf(1) = rf(isr) / beta
  rz(1) = rz(isr) / beta
  rzm1(1)  = 1.0d0 / rz(1)
  rf1_2(1) = sqrt(rf(1))
  rz1_2(1) = sqrt(rz(1))

  do i = isr+1,ier+1
     rf(i) = rf(i-1) * beta
     rz(i) = rz(i-1) * beta
     rzm1(i)  = 1.0d0 / rz(i)
     rf1_2(i) = sqrt(rf(i))
     rz1_2(i) = sqrt(rz(i))

! Calculate initial temperatures as well
     Tc(i) = T0*(rz(i)/AU)**(-p_T)

  enddo

  zgrid(isz) = 0.0d0
  do i 	= isz+1, iez
     zgrid(i)=zgrid(i-1) + dz
  enddo

  do i = isr, ier+1
     drfm1(i) = 1.0d0 / (rz(i)-rz(i-1))
  enddo

  do i = isr-1, ier
     drzm1(i) = 1.0d0 / (rf(i+1)-rf(i))
  enddo

! Set up sigma dot, due to winds and infall

  !call set_wind

  call set_accrete

  ! Set up surface density - iterates towards correct initial disk mass

  print*, 'Setting up disc surface density: Beginning iteration'

  mdisk = mdisk/solarmass
  rmax = rmax/AU

  if(sig_r==2) THEN
     sig_0_try = mdisk/(2.0*pi*log(rmax))
  ELSE
     sig_0_try = mdisk*(2.0-sig_r)/(2.0*pi*rmax**(2.0-sig_r))
  ENDif

  mdisktry = 0.0d0	

  sig_0_try = sig_0_try*solarmass/(AU*AU)

  tolerance = (mdisk-mdisktry)/mdisk
  do while(ABS(tolerance) > 0.001)

     mdisktry = 0.0d0

     do i = isr,ier 

  	if (rz(i).lt.rmax*AU) Then
           sigma(i) = sig_0_try*(rz(i)/AU)**(-sig_r)

  	Else
           sigma(i) = sig_0_try*(rz(i)/AU)**(-sig_r)
           sigma(i) = sigma(i)*DExp(-(rz(i)-rmax*AU)**2.0d0/(5.0d0*AU)**2.0d0)

  	Endif

	mdisktry = mdisktry + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma(i)

     enddo

     mdisktry = mdisktry/solarmass

     tolerance = (mdisk-mdisktry)/mdisk
     sig_0_try = sig_0_try*(1.0+ tolerance)

  enddo

  sigma_tot(:) = sigma(:)
  mdisk = mdisk*solarmass
  rmax = rmax*AU

  write(*,*) 'Disc has iterated initial mass of ',mdisktry

  close(12)

  ! Set surface density boundary conditions

  sigma(isr-1) = 0.0d0
  sigma(ier+1) = 0.0d0  

  !TODO - setup the layer properties here
  ! Calculate disc properties for this setup

  print*, 'Calculating other disc properties'
  call disc_properties
  write (*,*) '--- setup completed'


  return
end subroutine setup
