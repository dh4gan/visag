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

  tmig1 = mstar*mstar*aspectratio*aspectratio/(mdiscmig*mp(iplanet)*omegaK(iplanetrad(iplanet)))

  !print*, iplanetrad(iplanet),tmig1/yr, aspectratio,mp/mstar, mdiscmig/mstar, sigma(iplanetrad(iplanet))
  
end subroutine calc_typeI_migration
