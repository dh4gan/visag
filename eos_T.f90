      subroutine eos_T(rho,cs,T,kappa,mu,gamma,cp)
!-----------------------------------------------------------------------
! Reads and interpolates the equation of state tables to give
! temperatures, opacities ***USING T***
! Updated to fortran 90 PJC 22/05/2008
! Extended to opacities DHF 01/09/2009
!-----------------------------------------------------------------------

	use unitdata

      implicit none

      integer :: i,j
	real :: rho, cs,kappa,mu,gamma, T,cp
      real :: mcs,mcs1,mcs2,ccs,ccs1,ccs2,cs1,cs2
      real :: mmu,mmu1,mmu2,cmu,cmu1,cmu2,mu1,mu2
      real :: mkap,mkap1,mkap2,ckap,ckap1,ckap2,kap1,kap2 
      real :: mgamma,mgamma1,mgamma2,cgamma,cgamma1,cgamma2,gamma1,gamma2      
      

! Find the relevant records in the table...
! ... for rho
         if (rho < eostable(1,1,1)) rho = 1.01*eostable(1,1,1)
         i = 1
         do 
            if ((eostable(i,1,1) >= rho).or.(i == nrhopoints)) exit
            i = i + 1
         enddo

! ... and for internal energy
         if (T < eostable(1,1,2)) T = 1.01*eostable(1,1,2)
         j = 1
         do 
            if ((eostable(i-1,j,2) >= T).or.(j==nUpoints)) exit
            j = j + 1
         enddo
!         if(j==1) j=2
! Interpolate over the j value at i-1
         mcs1 = (eostable(i-1,j-1,3) - eostable(i-1,j,3))/ &
              (eostable(i-1,j-1,2) - eostable(i-1,j,2))
         ccs1 = eostable(i-1,j,3) - mcs1*eostable(i-1,j,2)
         cs1 = mcs1*T + ccs1

         mmu1 = (eostable(i-1,j-1,4) - eostable(i-1,j,4))/ &
              (eostable(i-1,j-1,2) - eostable(i-1,j,2))
         cmu1 = eostable(i-1,j,4) - mmu1*eostable(i-1,j,2)
         mu1 = mmu1*T + cmu1
		 
	 mkap1 = (eostable(i-1,j-1,6) - eostable(i-1,j,6))/ &
              (eostable(i-1,j-1,2) - eostable(i-1,j,2))
	 ckap1 = eostable(i-1,j,6) - mkap1*eostable(i-1,j,2)
	 kap1 = mkap1*T + ckap1
		 
	 mgamma1 = (eostable(i-1,j-1,5) - eostable(i-1,j,5))/ &
              (eostable(i-1,j-1,2) - eostable(i-1,j,2))
	cgamma1 = eostable(i-1,j,5) - mgamma1*eostable(i-1,j,2)
	gamma1 = mgamma1*T + cgamma1
				 
! Then interpolate over the j value at i
! Update j value as necessary
         j = 1
         do 
            if ((eostable(i,j,2) >= T).or.(j==nUpoints)) exit
            j = j + 1
         enddo

 !        if(j==1) j=2
         mcs2 = (eostable(i,j-1,3) - eostable(i,j,3))/ &
              (eostable(i,j-1,2) - eostable(i,j,2))
         ccs2 = eostable(i,j,3) - mcs2*eostable(i,j,2)
         cs2 = mcs2*T + ccs2
         
         mmu2 = (eostable(i,j-1,4) - eostable(i,j,4))/ &
              (eostable(i,j-1,2) - eostable(i,j,2))
         cmu2 = eostable(i,j,4) - mmu2*eostable(i,j,2)
         mu2 = mmu2*T + cmu2
		 
         mkap2 = (eostable(i,j-1,6) - eostable(i,j,6))/ &
              (eostable(i,j-1,2) - eostable(i,j,2))
         ckap2 = eostable(i,j,6) - mkap2*eostable(i,j,2)
         kap2 = mkap2*T + ckap2
		 
	 mgamma2 = (eostable(i,j-1,5) - eostable(i,j,5))/ &
              (eostable(i,j-1,2) - eostable(i,j,2))
         cgamma2 = eostable(i,j,5) - mgamma2*eostable(i,j,2)
         gamma2 = mgamma2*T + cgamma2

! Finally interpolate over i at the fractional j value
         mcs = (cs2 - cs1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ccs = cs2 - mcs*eostable(i,1,1)
         cs = mcs*rho + ccs

         mmu = (mu2 - mu1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cmu = mu2 - mmu*eostable(i,1,1)
         mu = mmu*rho + cmu
		 
	mkap = (kap2 - kap1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ckap = kap2 - mkap*eostable(i,1,1)
         kappa = abs(mkap*rho + ckap)
		 
	mgamma = (gamma2 - gamma1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cgamma = gamma2 - mgamma*eostable(i,1,1)
         gamma = mgamma*rho + cgamma

         cp = cs*cs/(T*gamma*(gamma-1))

      end subroutine eos_T
