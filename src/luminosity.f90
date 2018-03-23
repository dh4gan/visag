	subroutine luminosity(tot_lumin)

	use gravdata
	use unitdata
		
	implicit none

	integer :: i,j
	real(kind=8) :: tot_lumin,lumin,Blam(1000)  
	real(kind=8) :: lambda, dlambda, Teff
	
	do i = 1, 1000
  	spectrum(2,i) = 0.0d0
	Blam(i) = 0.0d0
	enddo
		  		  
	do i=isr, ier

	Teff = Tc(i)**4.0d0*tau(i)/(1.0d0+tau(i)**2.0d0)
  	Teff = Teff**0.25d0
			
    	lumin =  16.0d0/3.0d0*stefan*Tc(i)**4.0d0
    	lumin =  lumin*tau(i)/(1.0d0+tau(i)**2.0d0)
    	lumin =  lumin*(pi*rz(i+1)**2.0d0-pi*rz(i)**2.0d0)

    	tot_lumin = tot_lumin+lumin

    	lambda = 1.0d-4      ! 1 micron in cm
    	dlambda = 1.0d-4     ! 1 micron in cm
		    
    	do j = 1, 1000
      	Blam(j) = 2.0d0*planck*3.0d10**2.0d0/lambda**5.0d0
      	Blam(j) = Blam(j)/(DExp(planck*3.0d10/lambda/k_B/Teff)-1.0d0)
  
      	Blam(j) = Blam(j)*lambda*2.0d0*pi*rz(i)*(rz(i+1)-rz(i))

      	spectrum(1,j) = lambda

      	spectrum(2,j) = spectrum(2,j) + Blam(j)
      	spectrum(2,j) = spectrum(2,j)/4.0d0/pi/(140.0d0*pc)**2.0d0  ! distance of 140 pc

      	lambda = lambda + dlambda
    	enddo
			  
	enddo
			
	return
		      		   
	end subroutine luminosity
