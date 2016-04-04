subroutine midplane_properties(tauplus)
use unitdata, only:runmode
use gravdata, only:nmax
implicit none

real, dimension(nmax) :: tauplus

! Subroutine decides how to calculate midplane properties
! Given upper layer optical depth (tauplus)
! Depending on user selection

if(runmode=='g') then
    call midplane_properties_grav(tauplus)
else if (runmode=='Q') then
    call midplane_properties_grav_fixedQ(tauplus)
else
    call midplane_properties_fixedalpha(tauplus)
endif

end subroutine midplane_properties
