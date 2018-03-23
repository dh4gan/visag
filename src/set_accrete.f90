SUBROUTINE set_accrete
! Calculates sigma dot expected from mass infall rate
! (This version assumes a Gaussian distribution for sigma dot)

  use gravdata
  use unitdata, only: pi,yr,solarmass,AU
  use winddata
  
  implicit none

  integer :: i,in,out

  real :: inner,outer, term, integral
  real :: mean, sd

  ! Calculate discrete integral from adding together Gaussian over 3 * sd

  mean = 30.0*AU
  sd = 1.0*AU

  inner = mean - 3.0*sd
  outer = mean + 3.0*sd

  ! Find inner annulus

  i=1
  DO WHILE(rf(i) < inner)
     !print*, i,rf(i),inner
     i=i+1
  ENDDO

  in = i

  ! Now calculate integral
  DO WHILE(rf(i) < outer.and. i < ier-1)
  
     term= exp(-(rf(i)-mean)**2/(2.0*sd*sd))/sqrt(2.0*pi*sd*sd)
     integral = term*2.0*pi*rf(i)*(rf(i+1)-rf(i))
     i=i+1
  ENDDO
  out = i

!  print*, 'Integral is ', integral

! Now add integral to sigmadot

  DO i=in,out
     term =  exp(-(rf(i)-mean)**2/(2.0*sd*sd))/sqrt(2.0*pi*sd*sd)
     sigdot(i) = sigdot(i)+mdot_init*term/integral
  ENDDO


return
end subroutine set_accrete
