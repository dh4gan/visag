module unitdata
! Module stores all physical constants, EOS table,
! and other parameters relating to file storage

  integer,parameter :: itime = 80
  integer, parameter :: iprof = 83
  integer, parameter :: ispec = 78
  integer :: nrhopoints,nUpoints,counter
  real :: tdump, trun, tslice

  real(kind=8),parameter :: G = 6.67d-8
  real(kind=8), parameter :: AU = 1.496e13	 ! AU
  real(kind=8), parameter :: yr = 3.15e7	 ! year
  real(kind=8), parameter :: solarmass = 1.99e33 ! solar mass

  real(kind=8), parameter :: Qcrit = 1.5		 ! Critical Toomre Parameter	

  real(kind=8), parameter :: pi =3.1415926535
  real(kind=8), parameter :: stefan = 5.67d-5    ! Stefan-Boltzmann constant
  real(kind=8), parameter :: planck = 6.626d-27  ! Planck's constant
  real(kind=8), parameter :: pc = 3.08d18        ! parsec			
  real(kind=8), parameter :: k_B = 1.38d-16  	 ! Boltzmann constant
  real(kind=8), parameter :: m_H = 1.67d-24	 ! mass of hydrogen atom

  real(kind=8), allocatable, dimension(:,:,:) :: eostable

  character(len=100) :: fileno,prefix
end module unitdata
