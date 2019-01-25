subroutine calc_typeI_migration(iplanet, tmig1)
  !
  ! Computes the Type I migration timescale of a given planet
  !

  
  use planetdata
  use gravdata, only: mstar, omegaK,cs,sigma
  use unitdata, only: pi,yr

  implicit none

  integer, intent(in) :: iplanet
  real, intent(inout) :: tmig1
  real :: hp, aspectratio,mdiscmig


  ! Skip this calculation if planet not within disc
  if(iplanetrad(iplanet) < 0) then
     tmig1 = 1.0e30
  else
     
  
     hp = cs(iplanetrad(iplanet))/omegaK(iplanetrad(iplanet))
     aspectratio = hp/ap(iplanet)

     mdiscmig = pi*ap(iplanet)*ap(iplanet)*sigma(iplanetrad(iplanet))

     ! Baruteau et al (2013) expression
     tmig1 = mstar*mstar*aspectratio*aspectratio/(mdiscmig*mp(iplanet)*omegaK(iplanetrad(iplanet)))
  endif

end subroutine calc_typeI_migration
