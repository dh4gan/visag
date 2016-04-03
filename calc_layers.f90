
! ----------------------------------------------------------
! Subroutine to calculate t_cool
! ----------------------------------------------------------

subroutine calc_layers

  use gravdata
  use magdata
  use unitdata

  implicit none

  integer :: i

  IF(t/yr>t_activate) mag_switch = 1.0

  ! Begin by assuming no MRI activation
  ! Calculate disc properties assuming entirely self-gravitating

  tau_m(:) = 0.0
  sigma(:) = sigma_tot(:)

  CALL calc_grav(tau_m)

  IF(mag_switch==1.0) THEN
     ! Now look for MRI activation either fully or in the layer

     CALL calc_mag

     sigma_tot(:) = sigma(:) + sigma_m(:)

     ! Recalculate self-gravitating layer properties

     CALL calc_grav(tau_m)

     ! Otherwise, MRI switch not active and set layer's properties to zero
  ELSE
     nu_m(:) = 0.0
     sigma_m(:) = 0.0
     cs_m(:) = 0.0
     kappa_m(:) =0.0
     gamma_m(:) = 0.0
     mu_m(:) = 0.0
     tau_m(:) = 0.0
     sigma_tot(:) = sigma(:)
  ENDIF

  return


end subroutine calc_layers
