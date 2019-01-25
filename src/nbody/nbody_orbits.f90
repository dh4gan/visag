subroutine calc_orbit_from_vector
! Calculate the orbital parameters of all the bodies
! Also compute energy and angular momentum for error tracking

use embryodata
use eosdata, only: twopi,udist

implicit none

integer :: ibody,iembryo,ix
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

iembryo = ibody-1

if(embryo(iembryo)%finished==1) cycle

embryo(iembryo)%rmag = rmag(ibody)*udist

! Eccentricity first - calculate eccentricity (Laplace-Runge-Lenz) Vector
vdotr(:) = 0.0
do ix=1,3
   vdotr(:) = vdotr(:) + pos(ix,:)*vel(ix,:)
enddo

do ix=1,3
   eccvector(ix,:) = (vmag(:)*vmag(:)*pos(ix,:) -vdotr(:)*vel(ix,:))/gravparam - pos(ix,:)/rmag(:)
enddo

embryo(iembryo)%ecc = sqrt(eccvector(1,ibody)*eccvector(1,ibody) + &
    eccvector(2,ibody)*eccvector(2,ibody) + &
    eccvector(3,ibody)*eccvector(3,ibody))


! Semimajor axis
embryo(iembryo)%semimaj = angmag(ibody)*angmag(ibody)/&
    (gravparam*(1.0- embryo(iembryo)%ecc*embryo(iembryo)%ecc))

if(embryo(iembryo)%semimaj < 0.0) embryo(iembryo)%semimaj = abs(embryo(iembryo)%semimaj)

embryo(iembryo)%a = embryo(iembryo)%semimaj*udist

embryo(iembryo)%inc = 0.0

! Calculate the orbit's angles

   ! Inclination

   if(angmag(ibody)<small) cycle
   embryo(iembryo)%inc = acos(angmom(3,ibody)/ angmag(ibody))

   ! Longitude of the Ascending Node

   if (embryo(iembryo)%inc <small) then
      embryo(iembryo)%longascend = 0.0

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
         embryo(iembryo)%longascend = acos(nplane(1) / nmag);
      else
         embryo(iembryo)%longascend = 0.0
      endif

      if (nplane(2) < 0.0) embryo(iembryo)%longascend = twopi - embryo(iembryo)%longascend

   endif


   ! Calculate true anomaly

   !If orbit circular, no inclination, then use the position vector itself

   if (embryo(iembryo)%ecc < small .and. abs(embryo(iembryo)%inc) < small) then

      embryo(iembryo)%trueanom = acos(pos(1,ibody) / rmag(ibody));
      if (vel(1,ibody) < 0.0) embryo(iembryo)%trueanom = twopi - embryo(iembryo)%trueanom

      ! If orbit circular and inclination non-zero, then use the orbital plane vector
   else if (embryo(iembryo)%ecc < small) then

      ndotR = nplane(1)*pos(1,ibody) + nplane(2)*pos(2,ibody) + nplane(3)*pos(3,ibody)
      ndotR = ndotR / (rmag(ibody) * nmag);

      ndotV = nplane(1)*vel(1,ibody) + nplane(2)*vel(2,ibody) + nplane(3)*vel(3,ibody)

      embryo(iembryo)%trueanom = acos(ndotR);

      if (ndotV > 0.0) embryo(iembryo)%trueanom = twopi - embryo(iembryo)%trueanom

      ! For non-circular orbits use the eccentricity vector
   else

      edotR = eccvector(1,ibody)*pos(1,ibody) + eccvector(2,ibody)*pos(2,ibody) + eccvector(3,ibody)*pos(3,ibody)
      edotR = edotR / (rmag(ibody) * embryo(iembryo)%ecc);

      rdotV = vel(1,ibody)*pos(1,ibody) + vel(2,ibody)*pos(2,ibody) + vel(3,ibody)*pos(3,ibody)

     embryo(iembryo)%trueanom = acos(edotR);

      if (rdotV < 0.0)embryo(iembryo)%trueanom = twopi - embryo(iembryo)%trueanom
   endif

   ! Finally, calculate the longitude of periapsis - first calculate the argument of periapsis

   if (embryo(iembryo)%ecc > small) then

      edotn = eccvector(1,ibody)*nplane(1) + &
              eccvector(2,ibody)*nplane(2) + &
              eccvector(3,ibody)*nplane(3);

      edotn = edotn / (nmag * embryo(iembryo)%ecc);

      embryo(iembryo)%argper = acos(edotn);
      if (eccvector(3,ibody) < 0.0) embryo(iembryo)%argper = twopi - embryo(iembryo)%argper

   else

      embryo(iembryo)%argper = 0.0

   endif

enddo

end subroutine calc_orbit_from_vector

subroutine calc_vector_from_orbit
! Calculates body's position and velocity from orbital data

use embryodata
use eosdata, only: twopi, udist

implicit none

integer :: ibody,iembryo
real :: gravparam, semilatusrectum
!real :: a,e,i,long,om,nu
real,dimension(nbodies) :: rmag,vmag

real,dimension(nbodies) :: a,e,i,long,om,nu

a(:) = 0.0
e(:) = 0.0
i(:) = 0.0
long(:) = 0.0
om(:) = 0.0
nu(:) = 0.0

do ibody=2,nbodies
    iembryo=ibody-1

    a(ibody) = embryo(iembryo)%semimaj ! semimaj is in AU
    e(ibody) = embryo(iembryo)%ecc
    i(ibody) = embryo(iembryo)%inc
    long(ibody) = embryo(iembryo)%longascend
    om(ibody) = embryo(iembryo)%argper
    nu(ibody) = embryo(iembryo)%trueanom
    
    rmag(ibody) = embryo(iembryo)%semimaj * (1.0 - e(ibody) * e(ibody)) / (1.0 &
+ e(ibody) * cos(nu(ibody)))



! 2. Calculate position vector in orbital plane */

pos(1,ibody) = rmag(ibody)*cos(nu(ibody));
pos(2,ibody) = rmag(ibody) * sin(nu(ibody));
pos(3,ibody) = 0.0;

! 3. Calculate velocity vector in orbital plane */
semilatusrectum = abs(a(ibody) * (1.0 - e(ibody) * e(ibody)))
gravparam = totalmass

if (semilatusrectum > small) then

vmag(ibody) = sqrt(gravparam / semilatusrectum);

else

vmag(ibody) = 0.0;
endif


vel(1,ibody) = -vmag(ibody) * sin(nu(ibody));
vel(2,ibody) = vmag(ibody) * (cos(nu(ibody)) + e(ibody));
vel(3,ibody) = 0.0;

! 4. Begin rotations to correctly align the orbit
! Firstly, Rotation around z axis by -argument of Periapsis */

call rotate_Z(pos, nbodies, -1 * om);
call rotate_Z(vel, nbodies, -1 * om);

! Secondly, Rotate around x by -inclination */


call rotate_X(pos, nbodies, -1 * i);
call rotate_X(vel, nbodies, -1 * i);

! Lastly, Rotate around z by longitudeAscendingNode */

call rotate_Z(pos,nbodies,-1 * long);
call rotate_Z(vel,nbodies,-1 * long);

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





