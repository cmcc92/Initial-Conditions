
       program StellarWinds
       
       implicit double precision (a-z)

c---------------------------------------------------------------------
c                            Constants
c---------------------------------------------------------------------

       parameter (pi = 3.14159265359d0 ) 
       parameter (Rsun = 69.57d9) !cm
       parameter (CG = 6.67259d-8) !gravitational constant, cgs
       parameter (CMASS_EARTH = 5.972d27) !mass of earth, g
       parameter (CRADIUS_EARTH = 637.1d6) !radius of earth, cm
       parameter (CMP = 1.6726219d-24) !mass of proton, grams
       parameter (mmmars = 43.34d0) !molar mass of Mars (g/mol)

       parameter (CSOLARMASS_GRAM = 2.0d33) !solar mass to gram
       parameter (CYEAR_SEC = 31557600.0) !year to second
       parameter (CKM_CM = 1.0d05) !km to cm
       parameter (CCM_KM = 1.0d-5) !cm to km
       parameter (CAU_CM = 1.496d13) !AU to cm
       parameter (CG_KG = 0.001) !gram to kg
       
       parameter (co_2 = 44.0) !amu/q
       parameter (o_2 = 32.0) !amu/q
       parameter (nco2 = 4.0E5) !cm-3
       parameter (no2 = 3.6E6) !cm-3
       parameter (mco2 = 7.3E-23) !g
       parameter (mo2 = 5.3E-23) !g

c--------------------------------------------------------------------- 
c                          Common Blocks
c--------------------------------------------------------------------- 
       
       common/parkercomm/Rss,r2,omega,phi2
       common/planetBcomm/chapman2,rplanet2,rho2

c---------------------------------------------------------------------
c                        Reading Input File
c---------------------------------------------------------------------

       open(unit=11,file='example.in',status='old')
       rewind 11

       read (11,*) vr
       read (11,*) orbit_theta
       read (11,*) n0
       read (11,*) Bs0
       read (11,*) tempsw
       read (11,*) tempp
       read (11,*) chapman 
       read (11,*) rorbit 
       read (11,*) pmass
       read (11,*) pradius
       read (11,*) omega 
       read (11,*) Mdot 
       read (11,*) phi 
       read (11,*) rstellar 
       read (11,*) 

c---------------------------------------------------------------------
c                 Initial Conversion Calculations
c---------------------------------------------------------------------

       chapman2 = chapman*CKM_CM !cm (from km)
       rplanet = CRADIUS_EARTH*pradius !cm
       gravity = CG*pmass*CMASS_EARTH/(pradius*CRADIUS_EARTH)**2.0 !cm/s^2
       rorbit2 = rorbit*CAU_CM !cm
       rstar = rstellar*Rsun
       Rss = 2.6*Rstar !standard: 5 & avg. 2.6
       nsw = n0*(Rss/r)**2.0 !num density of stellar wind at planet 
       Mdot2 = Mdot*CSOLARMASS_GRAM/CYEAR_SEC !gram/s
       vr2 = vr*CKM_CM !cm/s
       phi2 = phi*pi/180.0
       omega2 = omega/(24*3600) !s^-1
       rhosw = nsw*CMP !mass density of stellar wind at planet
       rhop = nco2*mco2 + no2*mo2 !mass density of atmosphere
       mmp = CMP*C_ava !molar mass of proton
       mma = mmmars !molar mass of atmosphere
       
c---------------------------------------------------------------------
c                        Calling Subroutines
c---------------------------------------------------------------------
       call v_sound(temp_sw,molarP,cs)

       cs_sw = cs

       call v_sound(temp_p,molar_mars,cs)

       cs_p = cs
       
       call parker(orbit_theta,vrcgs,Bs0,Bx2,By2,Bz2,Btot)
       
       temp = Btot

       call planetBfield(vrcgs,temp,pBfield)

c---------------------------------------------------------------------
c                          Print Statements
c---------------------------------------------------------------------

       print*, cs_sw
       print*, cs_p
       print*, Bx2
       print*, By2
       print*, Bz2
       print*, Btot
       print*, pBfield

       stop
       end

c---------------------------------------------------------------------
c                           Subroutines
c---------------------------------------------------------------------
       
       subroutine v_sound(t,m,cs)
       implicit double precision (a-z)
       parameter (C_boltz = 1.3806d-16) !cgs
       parameter (C_ava = 6.022d23) !number/mol
       parameter (gama = 1.666666) 
       cs  = sqrt(gama*C_boltz*C_ava*t/m) 
       return
       end

c---------------------------------------------------------------------

       subroutine parker(theta,v,B,Bx2,By2,Bz2,Btot)
       implicit double precision (a-z)
       parameter (pi = 3.14159265359d0 ) 
       common/parkercomm/Rss,r2,omega2,phi2
       alpha = (90.0-theta)*pi/180.0
       theta2 = theta*pi/180.0 !convert to radians
       Br = B*(Rss/r2)**2.0
       Bphi = B*(Rss/r2)**2.0*(r2-Rss)*omega2*sin(theta2)/v
       Bx = Br*sin(theta2)*cos(phi)-(Bphi*sin(phi))
       By = Br*sin(theta2)*sin(phi)-(Bphi*cos(phi))
       Bz = Br*cos(theta2)
       alpha2 = -1.0*alpha
       Bx2 = cos(alpha2)*Bx-sin(alpha2)*Bz
       By2 = By
       Bz2 = sin(alpha2)*Bx2+cos(alpha2)*Bz
       Btot = sqrt(Bx2**2 + By2**2 + Bz2**2)
       return
       end

c---------------------------------------------------------------------

       subroutine planetBfield(v,B,pBfield)
       implicit double precision (a-z)
       parameter (pi = 3.14159265359d0 ) 
       common/planetBcomm/chapman2,rplanet2,rho2
       Bp = (8.0*pi*(B**2.0/(8.0*pi)+0.5*rho2*v**2.0))**0.5
       d = chapman2 + rplanet2
       pBfield = Bp*((rplanet2/d)**3.0)
       return
       end

c---------------------------------------------------------------------
c                              End
c---------------------------------------------------------------------
