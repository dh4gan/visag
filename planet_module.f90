module planetdata
! Contains all data relating to planets inside the disc

  ! Variables
  integer :: nplanet, nactive
  
  real(kind=8) :: rremove, p_create


  ! Arrays
  integer, allocatable,dimension(:)  :: alive_flag
  real, allocatable, dimension(:)  :: mp,ap, total_planet_torque, adot,tmig,tmigI
real, allocatable,dimension(:,:) :: lambdaI, lambdaII,fII
real, allocatable, dimension(:,:) :: torquei

end module planetdata
 
