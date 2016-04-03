      subroutine eos(rho,cs,temp,kappa,mu,gamma)
!-----------------------------------------------------------------------
! Reads and interpolates the equation of state tables to give
! temperatures, opacities ***USING c_s ***
! Updated to fortran 90 PJC 22/05/2008
! Extended to opacities DHF 01/09/2009
!-----------------------------------------------------------------------

     use unitdata

      implicit none

      integer :: i,j
      real :: rho, cs,temp,kappa,mu,gamma
      real :: mT,mT1,mT2,cT,cT1,cT2,T1,T2
      real :: mmu,mmu1,mmu2,cmu,cmu1,cmu2,mu1,mu2
      real :: mkap,mkap1,mkap2,ckap,ckap1,ckap2,kap1,kap2 
      real :: mgamma,mgamma1,mgamma2,cgamma,cgamma1,cgamma2,gamma1,gamma2
      

! Find the relevant records in the table...
! ... for rho

!      IF(t/yr> 999.0) print*, rho, eostable(1,1,1), cs, eostable(1,1,3)
         if (rho < eostable(1,1,1)) rho = 1.01*eostable(1,1,1)
         i = 1
         do 
            if ((eostable(i,1,1) >= rho).or.(i == nrhopoints)) exit
            i = i + 1
         enddo

! ... and for internal energy
         if (cs < eostable(1,1,3)) cs = 1.01*eostable(1,1,3)
         j = 1
         do 
            if ((eostable(i-1,j,3) >= cs).or.(j==nUpoints)) exit
            j = j + 1
         enddo
         IF(j==1) j=2

! Interpolate over the j value at i-1

         mT1 = (eostable(i-1,j-1,2) - eostable(i-1,j,2))/ &
              (eostable(i-1,j-1,3) - eostable(i-1,j,3))
         cT1 = eostable(i-1,j,2) - mT1*eostable(i-1,j,3)
         T1 = mT1*cs + cT1

         mmu1 = (eostable(i-1,j-1,4) - eostable(i-1,j,4))/ &
              (eostable(i-1,j-1,3) - eostable(i-1,j,3))
         cmu1 = eostable(i-1,j,4) - mmu1*eostable(i-1,j,3)
         mu1 = mmu1*cs + cmu1
		 
	 mkap1 = (eostable(i-1,j-1,6) - eostable(i-1,j,6))/ &
              (eostable(i-1,j-1,3) - eostable(i-1,j,3))
	 ckap1 = eostable(i-1,j,6) - mkap1*eostable(i-1,j,3)
	 kap1 = mkap1*cs + ckap1
		 
	 mgamma1 = (eostable(i-1,j-1,5) - eostable(i-1,j,5))/ &
              (eostable(i-1,j-1,3) - eostable(i-1,j,3))
	cgamma1 = eostable(i-1,j,5) - mgamma1*eostable(i-1,j,3)
	gamma1 = mgamma1*cs + cgamma1
				 
! Then interpolate over the j value at i
! Update j value as necessary
         j = 1
         do 
            if ((eostable(i,j,3) >= cs).or.(j==nUpoints)) exit
            j = j + 1
         enddo
         IF(j==1) j=2

         mT2 = (eostable(i,j-1,2) - eostable(i,j,2))/ &
              (eostable(i,j-1,3) - eostable(i,j,3))
         cT2 = eostable(i,j,2) - mT2*eostable(i,j,3)
         T2 = mT2*cs + cT2
         
         mmu2 = (eostable(i,j-1,4) - eostable(i,j,4))/ &
              (eostable(i,j-1,3) - eostable(i,j,3))
         cmu2 = eostable(i,j,4) - mmu2*eostable(i,j,3)
         mu2 = mmu2*cs + cmu2
		 
         mkap2 = (eostable(i,j-1,6) - eostable(i,j,6))/ &
              (eostable(i,j-1,3) - eostable(i,j,3))
         ckap2 = eostable(i,j,6) - mkap2*eostable(i,j,3)
         kap2 = mkap2*cs + ckap2
		 
	 mgamma2 = (eostable(i,j-1,5) - eostable(i,j,5))/ &
              (eostable(i,j-1,3) - eostable(i,j,3))
         cgamma2 = eostable(i,j,5) - mgamma2*eostable(i,j,3)
         gamma2 = mgamma2*cs + cgamma2

! Finally interpolate over i at the fractional j value
         mT = (T2 - T1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cT = T2 - mT*eostable(i,1,1)
         temp = mT*rho + cT

         mmu = (mu2 - mu1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cmu = mu2 - mmu*eostable(i,1,1)
         mu = mmu*rho + cmu
		 
	mkap = (kap2 - kap1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ckap = kap2 - mkap*eostable(i,1,1)
         kappa = mkap*rho + ckap
		 
	mgamma = (gamma2 - gamma1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cgamma = gamma2 - mgamma*eostable(i,1,1)
         gamma = mgamma*rho + cgamma

         ! Edit - remove opacity gap (Armitage et al 2001)
         kappa = 0.02*temp**0.8
      end subroutine eos
