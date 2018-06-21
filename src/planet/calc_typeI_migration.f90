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

  hp = cs(iplanetrad(iplanet))/omegaK(iplanetrad(iplanet))
  aspectratio = hp/ap(iplanet)

  mdiscmig = pi*ap(iplanet)*ap(iplanet)*sigma(iplanetrad(iplanet))

  ! Baruteau et al (2013) expression
  tmig1 = mstar*mstar*aspectratio*aspectratio/(mdiscmig*mp(iplanet)*omegaK(iplanetrad(iplanet)))

  ! Old Bate et al expression
  !tmig1 = aspectratio*mstar/(mp(iplanet)*omegaK(iplanetrad(iplanet)))
  

end subroutine calc_typeI_migration
