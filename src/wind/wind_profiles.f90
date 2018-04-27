! *******************************************
!
!   Font et al. and Alexander et al. wind profiles
!
! *******************************************

      Subroutine Wind_profiles(sigdot_diff,sigdot_dir, radius, Mstar, rdir_min)      

      real*8 sigdot_diff, sigdot_dir
      real*8 radius, rdir_min
      real*8 no, u1
      real*8 Rg, ng, cs
      real*8 mu, mH

      real*8 alphaB, C1, C2, A, B, D, Phi 
      real*8 Mstar, Msun
      real*8 G, pi

      pi = ACos(-1.0d0)
      G = 6.67d-8
      Msun = 1.989d33
      AU = 1.496d13

      mu = 1.35         ! ionised gas
      mH = 1.67d-24  

      cs = 1000000.0d0    ! 10 km/s sound speed in cm/s.
 
      C1 = 0.14
      C2 = 0.235 
      A = 0.3423
      B = -0.3612
      D = 0.2457

      alphaB = 2.6d-13
      Phi = 7.5d42

      If (Mstar.gt.0.0d0) Then

        Rg = G*Mstar/cs**2.0d0

        ng = C1*(3.0d0*Phi/(4.0d0*pi*alphaB*Rg**3.0d0))**0.5d0

        no = ng*(2.0d0/((radius/Rg)**(15.0d0/2.0d0)+(radius/Rg)**(25.0d0/2.0d0)))**(1.0d0/5.0d0)
   
        If ((radius.gt.0.1d0*Rg).and.(Rg/AU.gt.0.01d0)) Then
          u1 = cs*A*DExp(B*(radius/Rg - 0.1d0))*(radius/Rg - 0.1d0)**D
        Else
          u1 = 0.0d0
        EndIf 
  
        sigdot_diff = 2.0d0*no*u1*mu*mH
      Else
        sigdot_diff = 0.0d0
      EndIf   

      sigdot_dir = 2.0d0*C2*mu*mH*cs
      sigdot_dir = sigdot_dir*(Phi/(4.0d0*pi*alphaB*0.05*rdir_min**3.0d0))**0.5d0
      sigdot_dir = sigdot_dir*(radius/rdir_min)**(-2.42) 

      If (radius.lt.rdir_min) sigdot_dir = 0.0d0

      return
      end


