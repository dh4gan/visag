module planetdata
! Contains all data relating to planets inside the disc

  ! Variables
  integer, parameter :: npmax = 20
  integer :: np, nactive
  
  real(kind=8) :: rremove, p_create


  ! Arrays
  integer, allocatable,dimension(:) :: alive_flag
  real, dimension(npmax) :: mp,a

end module planetdata
 
