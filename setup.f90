!------------------------------------------------
! Subroutine to setup initial conditions and grid
!------------------------------------------------

subroutine setup

  use gravdata
  use planetdata
  use winddata
  use unitdata

  implicit none

  real(kind=8) :: rmax,sig_r
  real(kind=8) :: mdisk,mdisktry,dz
  real(kind=8) :: beta,rwindout,area, sig_0_try,tolerance

  integer :: i,ngrid


  open (unit=10,file='disc.par',form='formatted',status='old')

  read(10,*) prefix

  read (10,*) nrgrid
  read (10,*) nzgrid
  read(10,*) mstar
  read(10,*) mdisk
  read(10,*) sig_r
  read(10,*) T0
  read(10,*) p_T
  read (10,*) rin 
  read (10,*) rout
  read(10,*) rmax 
  read(10,*) zmax
  read (10,*) trun
  read (10,*) tdump	
  read (10,*) mdot_init
  read (10,*) rwind	
  read (10,*) mdot_wind	
  read (10,*) rremove

  prefix = TRIM(prefix)

  write(*,*) " "
  write(*,*) "-----------------------------------------------"
  write(*,*) "           SEMI ANALYTIC DISC CODE"
  write(*,*) "     Modified by D.Forgan, 29th April 2010       "
  write(*,*) "-----------------------------------------------"
  write(*,*) " "
  write(*,*) "-----------------------------------------------"
  write(*,*) " Input file: ./disc.par"
  write(*,*) "-----------------------------------------------"
  write(*,*) " "
  write(*,*) " "			
  write(*,*) "-----------------------------------------------"
  write (*,100) ' - nrgrid  = ',ngrid
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

  ! Convert variables to cgs units

  mdot_init = mdot_init * 6.2943d25
  mdot_wind = mdot_wind * 6.2943d25

  rin = rin*AU
  rout = rout*AU
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

  sigma(:) = 0.0
  Tc(:) = 0.0
  sigdot(:) = 0.0
! Set up source T according to input data
! (TODO)

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

  print*, 'Beginning sigma normalisation iteration'
  mdisk = mdisk/solarmass
  rmax = rmax/AU

  IF(sig_r==2) THEN
     sig_0_try = mdisk/(2.0*pi*log(rmax))
  ELSE
     sig_0_try = mdisk*(2.0-sig_r)/(2.0*pi*rmax**(2.0-sig_r))
  ENDIF

  mdisktry = 0.0d0	

  sig_0_try = sig_0_try*solarmass/(AU*AU)

  tolerance = (mdisk-mdisktry)/mdisk
  do while(ABS(tolerance) > 0.001)

     mdisktry = 0.0d0

     do i = isr,ier 

  	If (rz(i).lt.rmax*AU) Then
           sigma(i) = sig_0_try*(rz(i)/AU)**(-sig_r)

  	Else
           sigma(i) = sig_0_try*(rz(i)/AU)**(-sig_r)
           sigma(i) = sigma(i)*DExp(-(rz(i)-rmax*AU)**2.0d0/(5.0d0*AU)**2.0d0)

  	EndIf

	mdisktry = mdisktry + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma(i)

     enddo

     mdisktry = mdisktry/solarmass

     tolerance = (mdisk-mdisktry)/mdisk

     print*, mdisk,mdisktry,sig_0_try, tolerance
     sig_0_try = sig_0_try*(1.0+ tolerance)
     print*, sig_0_try

  enddo

  mdisk = mdisk*solarmass
  rmax = rmax*AU
  write(*,*) 'Disc has iterated initial mass of ',mdisktry


  close(12)

  ! Set surface density boundary conditions

  sigma(isr-1) = 0.0d0
  sigma(ier+1) = 0.0d0  

  write (*,*)
  write (*,*) '--- setup completed'

  snew(:) = 0.0
  call calc_grav(snew)

  return
end subroutine setup
