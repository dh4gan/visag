module nbodydata

  ! N Body variables


  integer, parameter :: inbodylog = 22
  real,parameter :: small = 1.0e-20
  real,parameter :: dampfac = 10.0
  real,parameter :: rsoft = 1.0e-5
  real,parameter :: tolerance = 1.0e-4

  real :: system_ang, system_energy, initial_system_ang,initial_system_energy
  real :: dE, dL, dt_nbody, totalmass, maxerror

  ! Body data

  integer :: nbodies
  real, allocatable, dimension(:,:) :: pos,vel,acc
  real, allocatable,dimension(:,:) :: newpos,newvel
  real,allocatable,dimension(:,:) :: angmom
  
  real,dimension(3) :: system_angmom,rcom,vcom,acom
  real,allocatable,dimension(:) :: mass, ekin,epot,etot,angmag
  real,allocatable,dimension(:) :: r, a, ecc, inc, longascend, argper, trueanom


contains

  subroutine get_magnitude(vector,magnitude)
    
    real,dimension(3),intent(in) :: vector
    real,intent(out):: magnitude
    
    magnitude = sqrt(vector(1)*vector(1)+vector(2)*vector(2)+vector(3)*vector(3))
    
  end subroutine get_magnitude

end module
  
