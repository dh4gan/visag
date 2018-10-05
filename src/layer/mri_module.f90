module magdata
  !---------------------------------------------------------------------------
  ! This module retains parameters and arrays used by the MRI turbulent layer
  ! (Feature in development!)
  !--------------------------------------------------------------------------
  
real, parameter :: tau_crit = 1.0
real, parameter :: alpha_m = 0.01
real, parameter :: Tcrit = 1000.0
real, parameter :: Toff =500.0
real, parameter :: t_activate = 5.0e4

real :: mag_switch
real, allocatable, dimension(:) :: sigma_m, nu_m,cs_m,kappa_m,sigma_tot
real, allocatable, dimension (:)  :: mu_m, gamma_m,tau_m, fullmri
end module magdata
