module unitdata
!-------------------------------------------------
! Module stores all physical constants, EOS table,
! and other parameters relating to file storage
!--------------------------------------------------
  
  integer,parameter :: itime = 80
  integer, parameter :: iprof = 83
  integer, parameter :: ispec = 78
  integer :: nrhopoints,nUpoints,snapshotcounter
  real :: tdump, trun, tslice,nfiles

  real(kind=8), parameter :: G = 6.67d-8           ! cgs
  real(kind=8), parameter :: AU = 1.496e13	   ! AU in cm
  real(kind=8), parameter :: yr = 3.15e7	   ! year in s
  real(kind=8), parameter :: solarmass = 1.99e33   ! solar mass
  real(kind=8), parameter :: mearth = 5.972e27     ! Earth mass
  real(kind=8), parameter :: mjup = 1.898e30       ! Jupiter mass
  real(kind=8), parameter :: msolyr = solarmass/yr

  real(kind=8), parameter :: Qcrit = 1.5	! Critical Toomre Parameter

  real(kind=8), parameter :: pi =3.1415926535
  real(kind=8), parameter :: twopi = 2.0*pi
  real(kind=8), parameter :: pibytwo = 0.5*pi
  real(kind=8), parameter :: stefan = 5.67d-5    ! Stefan-Boltzmann constant
  real(kind=8), parameter :: planck = 6.626d-27  ! Planck's constant
  real(kind=8), parameter :: pc = 3.08d18        ! parsec			
  real(kind=8), parameter :: k_B = 1.38d-16  	 ! Boltzmann constant
  real(kind=8), parameter :: m_H = 1.67d-24	 ! mass of hydrogen atom

  real(kind=8), allocatable, dimension(:,:,:) :: eostable

  character(len=1) :: runmode, layerchoice, tempchoice,nbodychoice,debug
  character(len=1) :: zerostring
character(len=6) :: snapshotformat
  character(len=100) :: fileno,prefix
  character(len=25), parameter :: paramfile = 'visag.params'


CONTAINS

  subroutine get_zero_padding_format(nfiles,zeroformat)
    implicit none

    integer,intent(in) :: nfiles
    integer :: nzeros
    character(1) :: zerostring
    character(6),intent(inout) :: zeroformat
    
    nzeros = int(log10(real(nfiles+1)))+2
    write(zerostring,'(I1)') nzeros
    zeroformat = "(I"//TRIM(zerostring)//"."//TRIM(zerostring)//")"
    
    return 
  end subroutine get_zero_padding_format
  
end module unitdata
