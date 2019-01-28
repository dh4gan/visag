module planetdata
! Contains all data relating to planets inside the disc

  real, parameter :: tdelay_planettorque = 1.0 ! Time delay for torques to activate in disc (years)

  ! Variables
  integer :: nplanet, nactive
  
  real(kind=8) :: rremove, p_create
character(100) :: planetfile
character(1) :: planetchoice

  ! Arrays
  integer, allocatable,dimension(:)  :: alive, iplanetrad, out_of_disc
  real, allocatable, dimension(:)  :: mp,ap
  real, allocatable, dimension(:) :: total_planet_torque, torque_term
  real, allocatable, dimension(:) :: adot,tmig,tmigI, fII
real, allocatable,dimension(:,:) :: lambdaI, lambdaII
real, allocatable, dimension(:,:) :: torquei

end module planetdata
 
