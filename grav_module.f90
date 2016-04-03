module gravdata
! This module contains all the properties of the self-gravitating disc layer
! Also contains mstar, t and dt for expediency

 integer :: nrgrid,nzgrid,isr,ier,isz,iez, nmax
 real :: rin,rout,zmax, p_T, T0
 real :: t,dt
 real(kind=8) :: mstar

  !	Zone centred arrays
  real, allocatable, dimension(:) :: rz,rzm1,rz1_2,drzm1
  !	Face centred arrays
  real, allocatable, dimension(:) :: rf,rf1_2,drfm1,zgrid

  ! Disc state variables (zone centred)
  real, allocatable, dimension(:) :: sigma, nu_tc, Tc, tau,omegaK,snew
  real, allocatable, dimension(:) :: cs, kappa, mu,gamma, tcool, alpha_g
  real,allocatable, dimension(:) :: T_source, sigdot,coolfunc,Q,Tnew,cp,heatfunc
  real, dimension(2,1000) :: spectrum

end module gravdata
