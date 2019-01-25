subroutine setup_planets
!
! Subroutine sets up planets to be added to disc
! Planet data read from separate input file
!
!

use gravdata
use planetdata
use unitdata
use nbodydata

implicit none

integer :: iplanet, ibody
character(6) :: filenumformat
character(100) :: outputfile

! Open planet file

open(10, file=planetfile,status='old')
! Read input data

read(10,*) nplanet

print*, 'There are ',nplanet, 'planets'
nactive = nplanet

allocate(mp(nplanet),ap(nplanet), ecc(nplanet),inc(nplanet))
allocate(longascend(nplanet), argper(nplanet), trueanom(nplanet))
allocate(alive(nplanet), iplanetrad(nplanet))
allocate(lambdaI(nplanet,nmax), lambdaII(nplanet,nmax))
allocate(fII(nplanet))
allocate(adot(nplanet),tmig(nplanet),tmigI(nplanet))
allocate(torquei(nplanet,nmax), torque_term(nmax), total_planet_torque(nmax))


alive(:) = 1
mp(:) = 0.0
ap(:) = 0.0
ecc(:) = 0.0
inc(:) = 0.0
longascend(:) = 0.0
argper(:) = 0.0
trueanom(:) = 0.0

iplanetrad(:) = 0

lambdaII(:,:) = 0.0
lambdaI(:,:) = 0.0
fII(:) = 0.0
adot(:) = 0.0
tmig(:) = 0.0
tmigI(:) = 0.0
torquei(:,:) = 0.0
total_planet_torque(:) = 0.0



do iplanet=1,nplanet
   read(10,*) mp(iplanet), ap(iplanet), ecc(iplanet), inc(iplanet), longascend(iplanet), argper(iplanet), trueanom(iplanet)
enddo
 

! Convert to correct units
mp(:) = mp(:)*mjup
ap(:) = ap(:)*AU


! If this is an N Body run, then create arrays for N body calculation
! Easier to do the N body calculation in separate arrays (which include star)
! N Body calculation units: (M=Msol, r=AU, t=1 yr/(2pi))

! Remember iembryo and ibody exclude/include star respectively

if(nbodychoice=='y') then
    nbodies = nplanet+1 ! Must include the star

    allocate(pos(3,nbodies),vel(3,nbodies),acc(3,nbodies), mass(nbodies))
    allocate(angmom(3,nbodies),angmag(nbodies))
    allocate(ekin(nbodies),epot(nbodies),etot(nbodies))
    allocate(newpos(3,nbodies),newvel(3,nbodies))

    pos(:,:) = 0.0
    vel(:,:) = 0.0
    acc(:,:) = 0.0

    angmom(:,:) = 0.0
    angmag(:) = 0.0
    ekin(:) = 0.0
    epot(:) = 0.0
    etot(:) = 0.0

    newpos(:,:) = 0.0
    newvel(:,:) = 0.0

    totalmass = totalmass/solarmass
    dt_nbody = 1.0e-3 ! Set arbitrary small timestep initially
    mass(1) = mstar/solarmass

    do ibody=2,nbodies
       iplanet = ibody-1
       mass(ibody) = mp(iplanet)/solarmass 
    enddo

    ! Debug lines - open N Body files to check orbits

    if(debug=='y') then
       
       ! Get character for body ID
       call get_zero_padding_format(nbodies,filenumformat)
       
       ! Open output files
       do ibody=2,nbodies          
          write(fileno, filenumformat) ibody          
          outputfile = TRIM(prefix)//".nbody."//TRIM(fileno)          
          open(ibody+inbodylog,file=outputfile, form="formatted")
       enddo
       
       ! Now set up log file
       outputfile = TRIM(prefix)//".nbody.log"

       open(inbodylog,file=outputfile,form="formatted")
    endif

    call calc_vector_from_orbit

endif


call find_planets_in_disc

do iplanet=1,nplanet
      print*, 'Planet ', iplanet, 'initially located at cell ', iplanetrad(iplanet)
      print*, 'Radius: ', rz(iplanetrad(iplanet))/AU
   enddo

end subroutine setup_planets
