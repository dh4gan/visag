! ********************************************
!   
!   X-ray wind from Owen
!
! ********************************************
       Subroutine xraywind(sigdotwind,xwind)
!
       real*8 sgdotwind, xwind
       real*8 a1, b1, c1, d1, e1, f1, g1
       real*8 temp1

       a1 = 0.15138
       b1 = -1.2182
       c1 = 3.4046
       d1 = -3.5717
       e1 = -0.32762
       f1 = 3.6064
       g1 = -2.4918

       If (xwind .lt. 0.7d0) Then
         sigdotwind = 0.0d0
       Else
         sigdotwind = a1*log10(xwind)**6.0d0 +b1*log10(xwind)**5.0d0
         sigdotwind = sigdotwind + c1*log10(xwind)**4.0d0    
         sigdotwind = sigdotwind + d1*log10(xwind)**3.0d0
         sigdotwind = sigdotwind + e1*log10(xwind)**2.0d0
         sigdotwind = sigdotwind + f1*log10(xwind) + g1
         sigdotwind = 10.0d0**sigdotwind

         temp1 = 6.0d0*a1*log(xwind)**5.0d0/(xwind**2.0d0*log(10.0d0)**7.0d0)
         temp1 = temp1 + 5.0d0*b1*log(xwind)**4.0d0/(xwind**2.0d0*log(10.0d0)**6.0d0)
         temp1 = temp1 + 4.0d0*c1*log(xwind)**3.0d0/(xwind**2.0d0*log(10.0d0)**5.0d0)
         temp1 = temp1 + 3.0d0*d1*log(xwind)**2.0d0/(xwind**2.0d0*log(10.0d0)**4.0d0)
         temp1 = temp1 + 2.0d0*e1*log(xwind)/(xwind**2.0d0*log(10.0d0)**3.0d0)
         temp1 = temp1 + f1/(xwind**2.0d0*log(10.0d0)**2.0d0)

         sigdotwind = sigdotwind*temp1*DExp(-1.0d0*(xwind/100.0d0)**10.0d0)
       EndIf

       Return
       End

