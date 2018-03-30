      subroutine eosread
!-----------------------------------------------------------------------
! Reads the equation of state tables
! Updated to fortran 90 PJC 22/05/2008
! Extended to opacities DHF 01/09/2009
!
! Note the columns (== k values) are as follows:
!   1) density
!   2) temperature
!   3) internal energy
!   4) mean molecular weight
!	5) Mass Weighted Opacity
!	6) Opacity
!
!  After recalculating the columns are:

! 1) density
! 2) temperature
! 3) sound speed
! 4) Mean Molecular Weight
! 5) Gamma
! 6) Opacity
!-----------------------------------------------------------------------

      use unitdata

      implicit none
      integer :: i,j,k,check

! Open eos data file and read in values

      open(50,file="myeos.dat", &
           status='old',iostat=check,action='read')
      if (check /= 0) then
         print*, 'Input EOS file not found'
         stop
      endif
      print*, "Reading equation of state table"
      
      read(50,*) nrhopoints, nUpoints

! Allocate input data array
      allocate(eostable(nrhopoints,nUpoints,6))
      do i = 1,nrhopoints
         do j = 1,nUpoints
            read(50,*) (eostable(i,j,k),k=1,6)
         enddo
      enddo
      close(50)

      print*, "Equation of state tables read in successfully"
      print*, 'nrhopoints, nUpoints, : ', nrhopoints, nUpoints	
      print*, "-----------------------------------------------"

	print*, "Calculating Gamma, sound speeds"

	do i =1,nrhopoints
	do j=1,nUpoints
	
!	Calculate gamma first (note column 3 still u)	
	eostable(i,j,5) = 1.0 + k_B*eostable(i,j,2)/(eostable(i,j,4)*m_H*eostable(i,j,3))
	 
!	Now calculate sound speed
	eostable(i,j,3) = eostable(i,j,5)*(eostable(i,j,5)-1)*eostable(i,j,3)
	

	enddo
	enddo	
			
	print*, 'Equation of state now fully initialised'

      return
      end subroutine eosread
        
