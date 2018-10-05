module gravdata
!---------------------------------------------------------
! This module contains all the properties of the disc
! Also contains mstar, t and dt for expediency
  !-------------------------------------------------------
  
 integer :: nrgrid,nzgrid,isr,ier,isz,iez, nmax
 real :: rin,rout,zmax, p_T, T_1AU
 real :: t,dt, alpha_visc, T_background
 real :: maxstep
 real(kind=8) :: mstar, mdot_init

  !	Zone centred arrays
  real, allocatable, dimension(:) :: rz,rzm1,rz1_2,drzm1
  !	Face centred arrays
  real, allocatable, dimension(:) :: rf,rf1_2,drfm1,zgrid

  ! Disc state variables (zone centred)
  real, allocatable, dimension(:) :: sigma, nu_tc, Tc, tau,omegaK,snew
  real, allocatable, dimension(:) :: cs, kappa, mu,gamma, tcool, alpha_g, H
  real,allocatable, dimension(:) :: T_source, sigdot,coolfunc,Q,Tnew,cp,heatfunc
  real, dimension(2,1000) :: spectrum

end module gravdata
