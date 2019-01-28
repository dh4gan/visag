subroutine calc_orbit_from_vector
! Calculate the orbital parameters of all the bodies
! Also compute energy and angular momentum for error tracking

  use nbodydata
  use planetdata
use unitdata, only: twopi

implicit none

integer :: ibody,iembryo,ix, iplanet
real :: gravparam
real :: rdotV, ndotR,ndotV,edotn,edotR, nmag
real,dimension(3) :: nplane
real,dimension(3,nbodies) :: eccvector
real,dimension(nbodies) :: vdotr,rmag,vmag

gravparam = totalmass

rmag(:) =sqrt(pos(1,:)*pos(1,:) + pos(2,:)*pos(2,:) + pos(3,:)*pos(3,:))
vmag(:) = sqrt(vel(1,:)*vel(1,:) + vel(2,:)*vel(2,:) + vel(3,:)*vel(3,:))

! Calculate orbital parameters - a,e,i

do ibody=2,nbodies

iplanet = ibody-1
  
if(alive(iplanet)==0) cycle

! Eccentricity first - calculate eccentricity (Laplace-Runge-Lenz) Vector
vdotr(:) = 0.0
do ix=1,3
   vdotr(:) = vdotr(:) + pos(ix,:)*vel(ix,:)
enddo

do ix=1,3
   eccvector(ix,:) = (vmag(:)*vmag(:)*pos(ix,:) -vdotr(:)*vel(ix,:))/gravparam - pos(ix,:)/rmag(:)
enddo

ecc(iplanet) = sqrt(eccvector(1,ibody)*eccvector(1,ibody) + &
    eccvector(2,ibody)*eccvector(2,ibody) + &
    eccvector(3,ibody)*eccvector(3,ibody))


! Semimajor axis
a(ibody) = angmag(ibody)*angmag(ibody)/&
    (gravparam*(1.0- ecc(iplanet)*ecc(iplanet)))

if(a(ibody) < 0.0) a(ibody) = abs(a(ibody))

! Calculate the orbit's angles

   ! Inclination

   if(angmag(ibody)<small) cycle
   inc(iplanet) = acos(angmom(3,ibody)/ angmag(ibody))

   ! Longitude of the Ascending Node

   if (inc(iplanet) <small) then
      longascend(iplanet) = 0.0

      nplane(1) = angmag(ibody)
      nplane(2) = 0.0
      nplane(3) = 0.0
      nmag = angmag(ibody)

   else

      nplane(1) = -angmom(2,ibody)
      nplane(2) = angmom(1,ibody);
      nplane(3) = 0.0;

      nmag = sqrt(nplane(1)*nplane(1) + nplane(2)*nplane(2) + nplane(3)*nplane(3));

      if(nmag>small) then
         longascend(iplanet) = acos(nplane(1) / nmag);
      else
         longascend(iplanet) = 0.0
      endif

      if (nplane(2) < 0.0) longascend(iplanet) = twopi - longascend(iplanet)

   endif


   ! Calculate true anomaly

   !If orbit circular, no inclination, then use the position vector itself

   if (ecc(iplanet) < small .and. abs(inc(iplanet)) < small) then

      trueanom(iplanet) = acos(pos(1,ibody) / rmag(ibody));
      if (vel(1,ibody) < 0.0) trueanom(iplanet) = twopi - trueanom(iplanet)

      ! If orbit circular and inclination non-zero, then use the orbital plane vector
   else if (ecc(iplanet) < small) then

      ndotR = nplane(1)*pos(1,ibody) + nplane(2)*pos(2,ibody) + nplane(3)*pos(3,ibody)
      ndotR = ndotR / (rmag(ibody) * nmag);

      ndotV = nplane(1)*vel(1,ibody) + nplane(2)*vel(2,ibody) + nplane(3)*vel(3,ibody)

      trueanom(iplanet) = acos(ndotR);

      if (ndotV > 0.0) trueanom(iplanet) = twopi - trueanom(iplanet)

      ! For non-circular orbits use the eccentricity vector
   else

      edotR = eccvector(1,ibody)*pos(1,ibody) + eccvector(2,ibody)*pos(2,ibody) + eccvector(3,ibody)*pos(3,ibody)
      edotR = edotR / (rmag(ibody) * ecc(iplanet));

      rdotV = vel(1,ibody)*pos(1,ibody) + vel(2,ibody)*pos(2,ibody) + vel(3,ibody)*pos(3,ibody)

     trueanom(iplanet) = acos(edotR);

      if (rdotV < 0.0)trueanom(iplanet) = twopi - trueanom(iplanet)
   endif

   ! Finally, calculate the longitude of periapsis - first calculate the argument of periapsis

   if (ecc(iplanet) > small) then

      edotn = eccvector(1,ibody)*nplane(1) + &
              eccvector(2,ibody)*nplane(2) + &
              eccvector(3,ibody)*nplane(3);

      edotn = edotn / (nmag * ecc(iplanet));

      argper(iplanet) = acos(edotn);
      if (eccvector(3,ibody) < 0.0) argper(iplanet) = twopi - argper(iplanet)

   else

      argper(iplanet) = 0.0

   endif

enddo

end subroutine calc_orbit_from_vector

subroutine calc_vector_from_orbit
! Calculates body's position and velocity from orbital data

use nbodydata
use planetdata
use unitdata, only: twopi, AU

implicit none

integer :: ibody,iplanet
real :: gravparam, semilatusrectum
real,dimension(nbodies) :: rmag,vmag

totalmass = sum(mass)

do ibody=2,nbodies
    iplanet =ibody-1
    
    rmag(ibody) = a(ibody) * (1.0 - ecc(ibody) * ecc(ibody)) / (1.0 &
+ ecc(ibody) * cos(trueanom(ibody)))


! 2. Calculate position vector in orbital plane */

pos(1,ibody) = rmag(ibody)*cos(trueanom(ibody));
pos(2,ibody) = rmag(ibody) * sin(trueanom(ibody));
pos(3,ibody) = 0.0;

print*, trueanom(ibody), rmag(ibody)
! 3. Calculate velocity vector in orbital plane */
semilatusrectum = abs(a(ibody) * (1.0 - ecc(ibody) * ecc(ibody)))
gravparam = totalmass

if (semilatusrectum > small) then

   vmag(ibody) = sqrt(gravparam / semilatusrectum);
   
else
   
   vmag(ibody) = 0.0;
endif


vel(1,ibody) = -vmag(ibody) * sin(trueanom(ibody));
vel(2,ibody) = vmag(ibody) * (cos(trueanom(ibody)) + ecc(ibody));
vel(3,ibody) = 0.0;

! 4. Begin rotations to correctly align the orbit
! Firstly, Rotation around z axis by -argument of Periapsis

call rotate_Z(pos, nbodies, -1 * argper)
call rotate_Z(vel, nbodies, -1 * argper)

! Secondly, Rotate around x by -inclination

call rotate_X(pos, nbodies, -1 * inc)
call rotate_X(vel, nbodies, -1 * inc)

! Lastly, Rotate around z by longitudeAscendingNode

call rotate_Z(pos,nbodies,-1 * longascend)
call rotate_Z(vel,nbodies,-1 * longascend)

enddo

end subroutine calc_vector_from_orbit

subroutine rotate_X(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in),dimension(nrows) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(:,:) = vector(:,:)
where(abs(angle)>1.0e-20)
newvector(1,:) = vector(1,:)
newvector(2,:) = vector(2,:)*cos(angle(:)) - vector(3,:)*sin(angle(:));
newvector(3,:) = vector(2,:)*sin(angle(:)) + vector(3,:)*cos(angle(:));
endwhere

vector(:,:) = newvector(:,:)

end subroutine rotate_X

subroutine rotate_Y(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in),dimension(nrows) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(:,:) = vector(:,:)
where(abs(angle)>1.0e-20)

newvector(1,:) = vector(1,:)*cos(angle(:)) + vector(3,:)*sin(angle(:));
newvector(2,:) = vector(2,:)
newvector(3,:) = -vector(1,:)*sin(angle(:)) + vector(3,:)*cos(angle(:));

endwhere

vector(:,:) = newvector(:,:)

end subroutine rotate_Y

subroutine rotate_Z(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in),dimension(nrows) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(:,:) = vector(:,:)

where(abs(angle)>1.0e-20)
newvector(1,:) = vector(1,:)*cos(angle(:)) - vector(2,:)*sin(angle(:));
newvector(2,:) = vector(1,:)*sin(angle(:)) + vector(2,:)*cos(angle(:));
newvector(3,:) = vector(3,:)
endwhere

vector(:,:) = newvector(:,:)

end subroutine rotate_Z





