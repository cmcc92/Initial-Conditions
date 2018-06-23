
       program StellarWinds
       
       implicit double precision (a-z)

c--------------------------------------------------------------------- 
c                          Common Blocks
c--------------------------------------------------------------------- 
       common/velocitycomm/gama,C_boltz,C_ava       
       common/parkercomm/pi,Rss,r2,omega,phi2

c---------------------------------------------------------------------
c                        Reading Input File
c---------------------------------------------------------------------

       open(unit=11,file='example.in',status='old')
       rewind 11

       read (11,*) num
       read (11,*) vr
       read (11,*) orbit_theta
       read (11,*) n0
       read (11,*) Bs0

c---------------------------------------------------------------------
c                            Constants
c---------------------------------------------------------------------

       pi = 3.14159265359d0
       Rsun = 69.57d9 !cm
       CG = 6.67259d-8 !gravitational constant, cgs
       CMASS_EARTH = 5.972d27 !mass of earth, g
       CRADIUS_EARTH = 637.1d6 !radius of earth, cm
       CMP = 1.6726219d-24 !mass of proton, grams
       rplanet = 9556.5d0 !km
       r = 1.0d0 !AU
       phi = 0.0d0 !degree
       Mdot = 4.0d-10 !solarmass/year
       temp_sw = 10.0d6 !temp of stellar wind, k 
       temp_p = 282.0d0 !temp of planet, K
       molar_mars = 43.34d0 !g/mol
       gama = 1.6666d0
       C_boltz = 1.3806d-16 !cgs
       C_ava = 6.022d23 !number/mol

c--------------------------------------------------------------------- 
c                       Conversion Constants
c--------------------------------------------------------------------- 

       CSOLARMASS_GRAM = 2.0d33 !solar mass to gram
       CYEAR_SEC = 365.25d0*24*3600 !year to second
       CKM_CM = 1000.0d0*100 !km to cm
       CAU_CM = 1.496d13 !AU to cm
       CCM_KM = 1.0d0/CKM_CM !cm to km
       CG_KG = 0.001 !gram to kg

c---------------------------------------------------------------------
c                           Ionosphere
c---------------------------------------------------------------------

       chap_layer = 61.2*100*1000 !cm (from km)
       co_2 = 44.0 !amu/q
       o_2 = 32.0 !amu/q
       nco2 = 4.0E5 !cm-3
       no2 = 3.6E6 !cm-3
       mco2 = 7.3E-23 !g
       mo2 = 5.3E-23 !g

c---------------------------------------------------------------------
c                 Initial Conversion Calculations
c---------------------------------------------------------------------

       rplanet2 = rplanet*CKM_CM !km to cm 
       gravity = CG*3.9*CMASS_EARTH/(1.5*CRADIUS_EARTH)**2.0 !cm/s^2
       r2 = r*CAU_CM !cm
       Rstar = 0.34*Rsun
       Rss = 2.6*Rstar !standard: 5 & avg. 2.6
       n_planet = n0*(Rss/r)**2.0 !num density of stellar wind at planet 
       Mdot2 = Mdot*CSOLARMASS_GRAM/CYEAR_SEC !gram/s
       vrcgs = vr*CKM_CM !cm/s
       phi2 = phi*pi/180.0
       omega = 0.44/(24*3600) !s^-1
       planet_rho = nco2*mco2 + no2*mo2
       molarP = CMP*C_ava
       molar_atmosphere = molar_mars
       
c---------------------------------------------------------------------
c                         Calling Functions
c---------------------------------------------------------------------
       call v_sound(temp_sw,molarP,cs)

       cs_sw = cs

       call v_sound(temp_p,molar_mars,cs)

       cs_p = cs
       
       call parker(orbit_theta,vrcgs,Bs0,Bx2,By2,Bz2,Btot)
       
c---------------------------------------------------------------------
c                          Print Statements
c---------------------------------------------------------------------

       print*, "Sound Speed of Atmosphere = ", cs_sw
       print*, "Sound Speed of Stellar Wind = ", cs_p
       print*, "Stellar Wind B-Field at Planet", Bx2
       print*, By2
       print*, Bz2
       print*, Btot

       stop
       end

c---------------------------------------------------------------------
c                             Functions
c---------------------------------------------------------------------
       
       subroutine v_sound(t,m,cs)
       implicit double precision (a-z)
       common/velocitycomm/gama,C_boltz,C_ava
       cs  = sqrt(gama*C_boltz*C_ava*t/m)
       return
       end

       subroutine parker(theta,v,B,Bx2,By2,Bz2,Btot)
       implicit double precision (a-z)
       common/parkercomm/pi,Rss,r2,omega,phi2
       alpha = (90.0-theta)*pi/180.0
       theta2 = theta*pi/180.0 !convert to radians
       Br = B*(Rss/r2)**2.0
       Bphi = B*(Rss/r2)**2.0*(r2-Rss)*omega*sin(theta2)/v
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









