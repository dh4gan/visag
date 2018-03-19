module planetdata
! Contains all data relating to planets inside the disc

  ! Variables
  integer :: np, nactive
  
  real(kind=8) :: rremove, p_create


  ! Arrays
  integer, allocatable,dimension(:)  :: alive_flag
     real, allocatable, dimension(:)  :: mplanet,aplanet, total_torque
     real, allocatable, dimension(:,:) :: torque

end module planetdata
 
