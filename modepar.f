ccccc6ccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
csoa     $ rhosv, rhob, catot, catote, sup, relh, frac, frabc, fraq, alpha)
      subroutine modepar(kcomp, ksol, imini, Nnat, rk, logsk, 
     $ rhosv, rhob, frombg, frbcbg, catot, catote, sup, relh, 
     $ frac, frabc, fraq, alpha)

c **********************************************************************************
c     Created by Alf Kirkevåg.
c **********************************************************************************

c     Prescribed dry lognormal number size distributions, accomodation coefficients
c     alpha, and hygroscopic swelling index (ksol=1 implies call of subroutine rhsub)  
c     for the clean background and for externally mixed SO4, BC and OC.
c     The index imin determines smallest radius used in the Mie calculations. 
c     We here also define the gridded concentrations of SO4+BC+OC for the tables, 
c     catot/catote, since this is different for each background mode kcomp, as well
c     as the gridded values for relh/sup, frac, frabc and fraq which are independent
c     on background mode number.   

c     Lag rent eksponensielle catot og catote som er "forutsigbare", s.a. en kan
c     forenkle deler av søkemekanismen i interpolasjonsrutinene i CAM-Oslo!?

      use commondefinitions

      implicit none

      INTEGER  imini, kcomp, ksol, isup
      REAL Nnat, rk, r0, rbcn, logsk, logs0, alpha
      REAL rhobc, rhooc, rhosv, rhob
      REAL catot(6), frac(6), frabc(6), fraq(6), relh(10), sup(9)
      REAL frombg(6), frbcbg(6)                                               !csoa
      REAL catote(16)

      Nnat=1.0    ! cm^(-3) normalized size distribution

      if(kcomp.ge.5) then  ! dummy array (defined but not used)
        catote=(/ 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
     $ 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10 /)
      endif
      if(kcomp.ge.1.and.kcomp.le.4) then  ! dummy array (defined but not used)
        catot=(/ 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10 /)
      endif

c     Median radius and log10(standard deviation) for the modes that are defined in CAM5-Oslo
c     (for those that do not exist there, the values given below counts)
      rk=originalNumberMedianRadius(kcomp)*1.e6
      logsk=log10(originalSigma(kcomp))
      if(kcomp.eq.1) then 
c        write(*,*) 'SO4(A/n), H2SO4 background for cond. of H2SO4'
cSOA        write(*,*) 'SO4(A/n), H2SO4 background for cond. of H2SO4 and SOA'
c       This is independent of iSOA (whether SOA is included or not)
        alpha=1.0
        ksol=1
        rhob=rhosv
        imini=1
        catote=(/ 1e-10, 1e-5, 2e-5, 4e-5, 8e-5, 1.5e-4, 3e-4,
     $  6e-4, 1.2e-3, 2.5e-3, 5e-3, 1e-2, 2e-2, 4e-2, 8e-2, 0.15 /)
      elseif(kcomp.eq.2) then 
c        write(*,*) 'BC(A/n), BC background for cond. of H2SO4'
        alpha=0.3
        ksol=1
        rhob=aerosol_type_density(2)
        imini=1   
        catote=(/ 1e-10, 1e-5, 2e-5, 4e-5, 8e-5, 1.5e-4, 3e-4,
     $  6e-4, 1.2e-3, 2.5e-3, 5e-3, 1e-2, 2e-2, 4e-2, 8e-2, 0.15 /)
      elseif(kcomp.eq.3) then     ! this mode is not defined/used in CAM-Oslo
c        write(*,*) 'OC(A/n), OC background for cond. of H2SO4' 
        alpha=0.7
        ksol=1
        rk=originalNumberMedianRadius(14)*1.e6
        logsk=log10(originalSigma(14))
        rhob=aerosol_type_density(3)
        imini=1   
        catote=(/ 1e-10, 1e-4, 2e-4, 4e-4, 8e-4, 1.5e-3, 3e-3,
     $  6e-3, 1.2e-2, 2.5e-2, 5e-2, 0.1, 2e-1, 0.4, 0.8, 1.5 /)
      elseif(kcomp.eq.4) then    
c        write(*,*) 'OC&BC(A/n), OC&BC background for cond. of H2SO4,' 
c        write(*,*) 'assuming OC is the basis for added BC and SO4' 
        alpha=0.5 ! between 0.3 for BC and 0.7 for OC
        ksol=1
        rhob=aerosol_type_density(3)
        imini=1   
        catote=(/ 1e-10, 0.01, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0,  ! std
     $   1.5, 2.5, 5., 10., 25., 50., 100., 500. /)*1.904e-3   ! std
      elseif(kcomp.eq.5) then    
c        write(*,*) 'SO4("Ait75"), H2SO4 background for cond/coag/Aq.'
        alpha=1.0
        ksol=1                            
        rhob=rhosv
        imini=13   
        catot=(/ 1.e-10, 5.e-4, 2.e-3, 0.01, 0.04, 0.15 /) 
      elseif(kcomp.eq.6) then    
c        write(*,*) '  MINACC, from AEROCOM'
        alpha=0.3
        ksol=1                            
        rhob=aerosol_type_density(4)
        imini=18      
        catot=(/ 1.e-10, 0.01, 0.05, 0.2, 0.8, 4.0 /)
      elseif(kcomp.eq.7) then    
c        write(*,*) '  MINCOA, from AEROCOM'
        alpha=0.3
        ksol=1                            
        rhob=aerosol_type_density(4)
        imini=20   
        catot=(/ 1.e-10, 0.02, 0.1, 0.5, 2.0, 8.0 /)
      elseif(kcomp.eq.8) then    
c        write(*,*) '  SSAIT, from AEROCOM'
        alpha=1.0
        ksol=1                            
        rhob=aerosol_type_density(5)
        imini=10      
cSS        catot=(/ 1e-10, 1e-4, 6e-4, 2.5e-3, 1e-2, 3.5e-2 /)
        catot=(/ 1.e-10, 5.e-4, 2.e-3, 0.01, 0.04, 0.15 /)  ! as for kcomp=5
      elseif(kcomp.eq.9) then    
c        write(*,*) '  SSACC, from AEROCOM'
        alpha=1.0
        ksol=1                            
        rhob=aerosol_type_density(5)
        imini=15     
cSS        catot=(/ 1.e-10, 0.005, 0.025, 0.1, 0.4, 2.0 /)
        catot=(/ 1.e-10, 0.01, 0.05, 0.2, 0.8, 4.0 /) ! as for kcomp=6
      elseif(kcomp.eq.10) then    
c        write(*,*) '  SSCOA, from AEROCOM'
        alpha=1.0
        ksol=1                            
        rhob=aerosol_type_density(5)
        imini=20
        catot=(/ 1.e-10, 0.02, 0.1, 0.5, 2.0, 8.0 /)
c     kcomp = 11, 12 and 14 are not used in CAM4-Oslo, just used for 
c     testing against kcomp = 1, 2 and 3 without condensate 
      elseif(kcomp.eq.0) then
c        write(*,*) '  soot (BC), fractal a-mode'
        alpha=0.3
        ksol=0
        rk=originalNumberMedianRadius(0)*1.e6
        logsk=log10(originalSigma(0))
        rhob=aerosol_type_density(2)    ! dummy, rhobcax used instead
        imini=1   
      else         
        write(*,*) 'modes 0 through 10 only'
        stop
      endif

c     define grid for tabulated optical parameters and CCN:
c     relative humidity, RH
      relh =(/ 0.0, 0.37, 0.47, 0.65, 0.75, 0.80, 
ct     $                    0.85, 0.90, 0.95, 0.98 /)
     $                    0.85, 0.90, 0.95, 0.995 /)
c     the fraction (internally mixed BC+OC)/(internally mixed BC+OC+SO4)
      frac=(/ 0.0, 0.1, 0.3, 0.5, 0.7, 0.999 /)
c     the fraction (internally mixed BC)/(internally mixed BC+OC)
      frabc=(/ 0.0, 0.01, 0.1, 0.3, 0.7, 0.999 /)
cM     the fraction (internally mixed SO4 from cloud processing)/(all internally mixed SO4)
c     the fraction (internally mixed (NH4)2SO4 from cloud processing)/(all internally mixed H2SO4 and (NH4)2SO4)
      fraq =(/ 0.0, 0.25, 0.50, 0.75, 0.85, 1.0   /)
c     supersaturation (for CCN calculations) in percents
      sup=(/ 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.8, 1.0 /)
c     supersaturation as a fraction
      do isup=1,9
        sup(isup)=1.0+0.01*sup(isup)
      enddo
csoa  the mass fraction OC/(OC + H2SO4) for the background aerosol of kcomp=1
      frombg=(/ 0.0, 0.1, 0.3, 0.5, 0.7, 0.999 /)
csoa  the mass fraction BC/(BC + OC) for the background aerosol of kcomp=4
      frbcbg=(/ 0.0, 0.1, 0.3, 0.5, 0.7, 0.999 /)
csoa

      return
      end 
