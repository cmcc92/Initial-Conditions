
       program Filecreater

c---------------------------------------------------------------------
c                            Constants
c---------------------------------------------------------------------

       implicit double precision (a-z) 
       parameter( pi = 3.14159265359d0 )
       parameter( Rsun = 69.57d9 ) !cm
       parameter( CG = 6.67259d-8 ) !gravitational constant, cgs
       parameter( CMASS_EARTH = 5.972d27) !mass of earth, g
       parameter( CRADIUS_EARTH = 637.1d6 ) !radius of earth, cm
       parameter( CMP = 1.6726219d-24 ) !mass of proton, grams
       parameter( rplanet = 9556.5d0 ) !km
       parameter( r = 1.0d0)  !AU
       parameter( phi = 0.0d0 ) !degree
       parameter( Mdot = 4.0d-10 ) !solarmass/year
       parameter( gama = 1.6666d0 ) !gamma
       parameter( temp_sw = 10.0d6 ) !temp of stellar wind, k 
       parameter( temp_p = 282.0d0 ) !temp of planet, K
       parameter( C_boltz = 1.3806d-16 ) !#cgs
       parameter( C_ava = 6.022d23 ) !number/mol
       parameter( molar_mars = 43.34d0 ) !g/mol

c--------------------------------------------------------------------- 
c                       Conversion Constants
c--------------------------------------------------------------------- 

       parameter( CSOLARMASS_GRAM = 2.0d33 ) !solar mass to gram
       parameter( CYEAR_SEC = 365.25d0*24*3600 ) !year to second
       parameter( CKM_CM = 1000.0d0*100 ) !km to cm
       parameter( CAU_CM = 1.496d13 ) !AU to cm
       parameter( CCM_KM = 1.0d0/CKM_CM ) !cm to km
       parameter( CG_KG = 0.001 ) !gram to kg

c---------------------------------------------------------------------
c                           Ionosphere
c---------------------------------------------------------------------

       parameter( chap_layer = 61.2*100*1000 ) !cm (from km)
       parameter( co_2 = 44.0 ) !amu/q
       parameter( o_2 = 32.0 ) !amu/q
       parameter( nco2 = 4.0E5 ) !cm-3
       parameter( no2 = 3.6E6 ) !cm-3
       parameter( mco2 = 7.3E-23 ) !g
       parameter( mo2 = 5.3E-23 ) !g

c--------------------------------------------------------------------- 
c                          Common Blocks
c--------------------------------------------------------------------- 

      common/constants/pi,Rsun

c---------------------------------------------------------------------
c                        Reading Input File
c---------------------------------------------------------------------

       open(unit=11,file='initialinput2',status='old')
       rewind 11

       read (11,*) num
       read (11,*) vr
       read (11,*) theta
       read (11,*) n0
       read (11,*) Bs0

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
       Omega = 0.44/(24*3600) !s^-1
       rho = nco2*mco2 + no2*mo2
       molarP = CMP*C_ava
       molar_atmosphere = molar_mars
       
       call v_sound(temp_sw,molar_mars,gama,C_boltz,C_ava,cs) 
       
       cs_sw = cs

       call v_sound(temp_p,molar_mars,gama,C_boltz,C_ava,cs)

       cs_p = cs

c---------------------------------------------------------------------
c                        Print Statements
c---------------------------------------------------------------------

       print*, num
       print*, vr
       print*, theta
       print*, n0
       print*, Bs0
       print*, cs_sw
       print*, cs_p

       stop
       end

c---------------------------------------------------------------------
c                       Functions
c---------------------------------------------------------------------
       
       subroutine v_sound(T,M,G,boltz,ava,cs)
       double precision T,M,G,boltz,ava,cs
       cs = sqrt(G*boltz*ava*T/M)
       return
       end
