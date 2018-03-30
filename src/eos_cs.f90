      subroutine eos_cs(rho,cs,temp,kappa,mu,gamma,cp)
!-----------------------------------------------------------------------
! Reads and interpolates the equation of state tables to give
! temperatures, opacities ***USING sound speed***
! Updated to fortran 90 PJC 22/05/2008
! Extended to opacities DHF 01/09/2009
!-----------------------------------------------------------------------

	use unitdata

      implicit none

      integer :: i,j
	real, intent(inout) :: rho, cs,kappa,mu,gamma, temp,cp
      real :: mtemp,mtemp1,mtemp2,ctemp,ctemp1,ctemp2,temp1,temp2
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

! ... and for sound speed         

         if (cs < eostable(1,1,3)) cs = 1.01*eostable(1,1,3)
         j = 1
         do 
            if ((eostable(i-1,j,3) >= cs).or.(j==nUpoints)) exit
            j = j + 1
         enddo

!         if(j==1) j=2
! Interpolate over the j value at i-1
         mtemp1 = (eostable(i-1,j-1,2) - eostable(i-1,j,2))/ &
              (eostable(i-1,j-1,3) - eostable(i-1,j,3))
         ctemp1 = eostable(i-1,j,2) - mtemp1*eostable(i-1,j,3)
         temp1 = mtemp1*cs + ctemp1

        ! print*, cs, j, mtemp1,ctemp1,temp1

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
				 
! then interpolate over the j value at i
! Update j value as necessary
         j = 1
         do 
            if ((eostable(i,j,3) >= cs).or.(j==nUpoints)) exit
            j = j + 1
         enddo

 !        if(j==1) j=2
         mtemp2 = (eostable(i,j-1,2) - eostable(i,j,2))/ &
              (eostable(i,j-1,3) - eostable(i,j,3))
         ctemp2 = eostable(i,j,2) - mtemp2*eostable(i,j,3)
         temp2 = mtemp2*cs + ctemp2
         
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
         mtemp = (temp2 - temp1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ctemp = temp2 - mtemp*eostable(i,1,1)
         temp = mtemp*rho + ctemp

         mmu = (mu2 - mu1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cmu = mu2 - mmu*eostable(i,1,1)
         mu = mmu*rho + cmu
		 
	mkap = (kap2 - kap1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ckap = kap2 - mkap*eostable(i,1,1)
         kappa = abs(mkap*rho + ckap)
		 
	mgamma = (gamma2 - gamma1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cgamma = gamma2 - mgamma*eostable(i,1,1)
         gamma = mgamma*rho + cgamma

         cp = cs*cs/(temp*gamma*(gamma-1))

      end subroutine eos_cs
