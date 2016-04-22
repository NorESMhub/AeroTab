
      program AeroTab   ! Program for making Aerosol look-up tables for CAM5-Oslo

c **********************************************************************************
c     Created by Alf Kirkevåg. The code is originally based on the method developed 
c     and described by Kirkevåg, A., Iversen, T., and Dahlback, A. (1999): On radiative 
c     effects of black carbon and sulphate aerosols. Atmos. Environ. 33, 2621-2635. 
c **********************************************************************************

c     This program defines initial size distributions (at the point of emission) 
c     and microphysical properties, such as hygroscopicity and wevelength dependent 
c     refractive indices, then calculates the modified aerosol size distributions 
c     (after aerosol processing), and finally it performs Mie and CCN or dry size 
c     parameter calculations for the clean or internally mixed aerosol modes as a 
c     function of relative humidity (or supersaturation for CCN) and added internally 
c     mixed BC, OC and sulfate. The output is a range of look-up tables (*.out) for 
c     use in CAM5-Oslo (with only some CAM-Oslo/CAM4-Oslo functionality retained). 
c
c     References for CAM-Oslo (based on CAM3):
c     Kirkevåg, A., Iversen, T., Seland, Ø., Debernard, J.B., Storelvmo, T., and 
c      Kristjánsson, J.E. (2008) Aerosol-cloud-climate interactions in the climate 
c      model CAM-Oslo. Tellus, 60A, 492-512.
c     Seland, Ø., T. Iversen, A. Kirkevåg, and T. Storelvmo (2008) Aerosol-climate 
c      interactions in the CAM-Oslo atmospheric GCM and investigations of associated 
c      shortcomings. Tellus, 60A, 459-491.
c     Reference for CAM4-Oslo (based on CAM4):
c     Kirkevåg, A., T. Iversen, Ø. Seland, C. Hoose, J. E. Kristjánsson, H. Struthers, 
c      A. Ekman, S. Ghan, J. Griesfeller, D. Nilsson, and M. Schulz: Aerosol-climate 
c      interactions in the Norwegian Earth System Model - NorESM1-M, Geosci. Model Dev., 
c      6, 207-244, doi:10.5194/gmd-6-207-2013, 2013. 
c     CAM5-Oslo is at present under development. Reference for nucleation and simplest
c     SOA treatment (assuming SOA->SO4): Makkonen, R., Seland, Ø., Kirkevåg, A., Iversen, 
c     T., and Kristjánsson, J. E.: Evaluation of aerosol number concentrations in NorESM 
c     with improved nucleation parameterisation, Atmos. Chem. Phys. 14, 5127-5152, 
c     doi:10.5194/acp-14-5127-2014, 2014. 

C ===================================================================================
C Notes on past and present code development:
C     Since the CAM4-Oslo look-up tables were made, some parts of the code have been
C     slightly modified: The diffusion coefficient and mean free path for sulfuric acid 
C     (H2SO4) have been updated, see constize.f. Therefore, do not expect to reproduce 
C     the old look-up tables exactly as they were (small changes).  
C
C     One inconsistency compared with the CAM4-Oslo life cycle model description in 
C     Kirkevåg et al. (2013), is that the OC(n) mode (kcomp=3 without condensed SO4), 
C     is still used when excessive OM mass (exceeding the max table values defined in
C     modepar.f) is lumped in the model. Should we remove this mode and rather lump 
C     mass to a clean OM&BC(Ait) mode (kcomp=4 without condensed SO4)? (See below).
cSOA
c     April 2013:
C     Due to new SOA treatment, the OC(n) mode may be needed anyway. Possible changes 
C     to accomodate this may be: add condensed SOA to BC(n/Ait), OC(n/Ait), and (as 
C     proposed by R. Makkonen) SO4(n/Ait). So far only the SO4(n/Ait) mode has been 
C     included in this code. 
cSOA   
C     A possible simplification to think about for CAM5-Oslo: 
C     let DU(c) and SS(c) be externally mixed only. 
C     Additional functionality:
C     Consider to include for ib=31: abs550 for each component (in aerocomk*.out)!
cM
C     August 2014: Removing dependencies on molar weights (Ms, Msv and Mso4) to 
C     facilitate corresponding simplifications in CAM5-Oslo. 
cM
c13/cn 
C     April 2015: remove modes (kcomp=) 11,12,14, and rename mode 13 to 0.
C     r12 is renamed rbcn.
c13/cn
csoa
c     July 2015:
c     Including new SOA treatment, allowing for condensation of VOC/SOA onto all 
c     background modes kcomp= 1 - 10, but treated as OM coagulate for modes 5-10
c     (the same way that H2SO4 condensate is treated as coagulate for these modes).
c     We still keep kcomp=3, since it may be used in parts of the code (for lumping
c     of overshooting mass w.r.t. upper ceiling in the look-up tables), and since
c     there is still no need for a new mode to fill its "place". 
c     New treatment for kcomp=1 & 4: fombg and fbcbg are now mass fractions of OM 
c     or BC in the background aerosol (not radius dependent), instead of using the
c     trick for kcomp=4 of redifining fac. fac has now the same meaning for all 
c     modes. 
c     
csoa
c     NB: det ser ut til at hygroskopisk vokst av OM&BC(Ait) før bare var bestemmes           !!!
c     av OM, ikke BC! Rettet nå!                                                              !!! 
c     Forslag til ny feature (relatert til det ovenfor): anta coating mhp RH-vokst            !!!
c     på samme moder som dette skjer med hensyn på CCN aktivering i CAM5-Oslo!                !!!
C ===================================================================================

      implicit none

ccccc6ccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      INTEGER  kcom, iband
      INTEGER  i, imax, imini, imaxi, is, icl, ib, ibcam, ictot, kcomp,
     $ kc, ksol, itot, ifbc, ifac, ifaq, irelh, iopt, isup, ismolar, 
     $ ismolarh, irh, itilp, 
     $ ifombg, ifbcbg                                                  ! csoa
      INTEGER isup1, isup2, irelh1, irelh2, ictot1, ictot2, ifac1, 
     $ ifac2, ifbc1, ifbc2, ifaq1, ifaq2, ictote, ictote1, ictote2,
     $ ifombg1, ifombg2, ifbcbg1, ifbcbg2                              ! csoa
      REAL r(0:100), rp(0:100), fki(-1:100), fracdim(0:100),  
     $ dndlrk0(0:100), dndlrk(0:100), dndlrkny(0:100), 
     $ vbci(0:100), voci(0:100), vsi(0:100), vai(0:100), vssol(0:100), 
     $ vbcsol(0:100), vocsol(0:100), vasol(0:100), vw(0:100) 
      REAL rk, r0, rbcn, rcoag, d, ntot, Nnatk, Nnat, fcondk, 
     $ fcoagk, faqk, logsk, logs0, 
cM     $ Ms, Mso4, Msv, 
     $ rhos, rhosv, rhoc2, rhobc, rhooc, rhob, rhow, th, mfv, diff,  
     $ Cac, Cabc, Caoc, Cas1, Cas2, Cas3, Caso4, Ctot, Cdry, dCtot, cat, 
cM     $ fr, fac, fabc, faq, rh, supers, rsup, numb, bcint, CCN, cintbg, 
     $ fac, fabc, faq, rh, supers, rsup, numb, bcint, CCN, cintbg, 
     $ cintsu, cintsc, cintsa, cintbc, cintoc, cintbg05, cintsu05, 
     $ cintsc05, cintsa05, cintbc05, cintoc05, cintbg125, cintsu125, 
     $ cintsc125, cintsa125, cintbc125, cintoc125, aaero, aaeros, 
     $ aaerol, vaero, vaeros, vaerol, bclt05, bcgt125, lambda, alpha, 
     $ fombg, vombg, fbcbg, vbcbg, eps                                 ! csoa
      REAL frombg(6), frbcbg(6)                                        ! csoa
      REAL catot(6), frac(6), frabc(6), fraq(6), relh(10), sup(9)
      REAL catote(16)
      REAL omega(31), gass(31), bext(31), babs(31), kext(31)
      REAL xlam(31), xlami(32), xlamb(31), xlame(31),
     $ fband(31), fb(16)
      REAL Ctot0, Dm(0:100), Dmp(0:100), K12(0:101), Kp12(0:101),
     $ K12oc(0:101), Kp12oc(0:101), K12so4(0:101), Kp12so4(0:101), 
     $ Ctotnull                               
      REAL rcoag_so4n, rcoag_bcn, rcoag_ocn
      REAL xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss
cSOA  Added May & December 2013
      REAL diffsoa, thsoa, mfvsoa, Dmsoa(0:100), Dmpsoa(0:100), casoa
csoa      INTEGER iSOA
cSOA
      REAL pi, e, testnumb 
      COMPLEX cref(5,31)
      PARAMETER  (pi=3.141592654, e=2.718281828)
      PARAMETER (eps=1.e-50)

c     Assumed radius for coagulating (fine mode) particles:  
c     (a common coagulation radius rcoag=0.04 was used originally)
      PARAMETER (rcoag_so4n = 0.0118)  ! rk for kcomp=1
      PARAMETER (rcoag_bcn  = 0.0118)  ! rk for kcomp=2
      PARAMETER (rcoag_ocn  = 0.04)    ! rk for kcomp=3

c     Do not modify the following input:
c     number of iterations in the Smolarkiewicz advection scheme, ismolar, and
c     a different ismolarh for hygroscopic growth than for dry distributions:
      PARAMETER  (ismolar=2, ismolarh=3)
c     no aportioning between modes (i.e., all material is added onto the same mode): 
c     (this is relevant to modify only when this code is part of a larger multimodal
c     scheme, i.e. for distribution of internally mixed mass onto more than one mode
c     at the time).
      PARAMETER  (Nnatk=1.0, fcondk=1.0, fcoagk=1.0, faqk=1.0)

c     Modify the following input to create different sets of look-up tables:
c     Let iopt=1 for optics tables, iopt=0 for CCN (CAM-Oslo with diagnostic CDNC)
c     --> ccnk*.out, or size distribution calculations (CAM4-Oslo and CAM5-Oslo 
c     with the prognostic CDNC scheme):
cNB   the combination iopt=0 and itilp=0 (giving CCN as function of S) is not updated: 
cNB   do not use without cleaning up and checking first! (something wrong for kcomp=8-10,
cNB   at least). Look-up tables for this combination is not needed in any CAM4-Oslo or 
cNB   CAM5-Oslo version... 
      iopt=1
c     Lognormal mode fitting (itilp=1) or not (itilp=0) (requires iopt=0)
c     --> logntilp*.out (and nkcomp.out for dry size distributions). To ensure
c     that only logntilp*.out and not ccnk*.out tables are produced, we let:
      itilp=1-iopt  
c     Outout for iopt=1 --> lwkcomp*.out or kcomp*.out, aerodryk*.out, 
c     aerocomk*.out, and nkcomp*.out (for size distributions for all RH).
c     SW: ib=29 (ave.=>12) SW "bands" (CAMRT), or
c     SW: ib=31 (ave.=>14) (RRTMG) (Added November 2013), or
c     LW: ib=19 (ave.=>16) (RRTMG) (Added November 2013):
      ib=31
cSOA  Added December 2013
c     SOA may be internally mixed with the SO4(ait) mode (1) or not (0). 
c     iSOA=0 in CAM4-Oslo/NorESM1 (e.g., Kirkevåg et al., 2013)
csoa      iSOA=1
cSOA

C     Initialization and calculations of look-up tables starts here...       

      if(ib.eq.29) then
          write(*,*) 
     $ 'Note: for aerocomk*.out, aerodryk*.out or SW RRTMG, use ib=31'    
      endif

c     Define spectral bands and spectral solar fluxes (at TOA) to be used 
c     in Chandrasekhar averaging of the optical parameters (in sizemie)
      call specbands(ib, xlami, xlam, xlamb, xlame, fband, fb, ibcam)

c     Define constants and parameters for calculations of size distributions   ! Move some of this to modepar.f !!!
      call constsize(d, imax, imaxi, r, rp, r0, rbcn, logs0,
     $ rhobc, rhooc, rhos, rhosv, rhoc2, rhow, 
     $ bcint, fracdim, diff, th, mfv, diffsoa, thsoa, mfvsoa)

c     The main loop over aerosol mode number for background modes, kcomp=1,10, 
c     plus kcomp=0, for the fractal BC(ac) mode (use itot=0). For this mode no 
c     lognormal mode fitting is needed (it is assumed to be hydrophobic and 
c     therefore not giving any CCN or CDNC contribution).  

      do kcomp=0,10   ! for look-up tables, kcomp=0,10 (only kcomp=1-10 needed for logntilp*.out)

       if(kcomp.eq.0) then
          itot = 0     ! not subject to added mass by condensation etc.
       else
          itot = 1     ! subject to added mass by condensation etc.
       endif

c     Set parameters for prescribed initial dry lognormal size 
c     distributions, and grid for tabulated optical parameters (or CCN)
ccccc6ccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
csoa      call modepar(kcomp, ksol, imini, Nnat, rk, logsk, 
csoa     $ rhosv, rhob, catot, catote, sup, relh, frac, frabc, fraq, alpha)
      call modepar(kcomp, ksol, imini, Nnat, rk, logsk, rhosv, rhob, 
     $ frombg,frbcbg,catot, catote, sup, relh, frac, frabc, fraq, alpha)

c     drydist calculates dry background mode size distribution, dndlrk0,  
c     and dry aerosol contribution to mass concentration, Ctot0 (ug/m**3). 
      call drydist(kcomp, Nnat, imini, imax, d, r, rk,  
     $ logsk, logs0, rhob, bcint, pi, dndlrk0, ntot, Ctot0)

      Ctotnull=Ctot0

c     Diffusion and coagulation coeffecients for the respective background
c     mode are then calculated. We here assume that all n/Aitken-modes that 
c     coagulate on a mineral or sea-salt background all are monodisperse.
c
c     Diffusion coefficients Dm for H2SO4
      call condsub (r, imax, diff, mfv, th, alpha, Dm)     
      call condsub (rp, imax, diff, mfv, th, alpha, Dmp)     
c     Diffusion coefficients Dm for SOA
      call condsub (r, imax, diffsoa, mfvsoa, thsoa, alpha, Dmsoa)     
      call condsub (rp, imax, diffsoa, mfvsoa, thsoa, alpha, Dmpsoa)     
c     Coagulation coefficients K12 for H2SO4
      call coagsub (r, imax, rcoag_so4n, rhob, rhosv, K12so4)   
      call coagsub (rp, imax, rcoag_so4n, rhob, rhosv, Kp12so4)   
c     Coagulation coefficients K12 for BC
      call coagsub (r, imax, rcoag_bcn, rhob, rhobc, K12)   
      call coagsub (rp, imax, rcoag_bcn, rhob, rhobc, Kp12)   
c     Coagulation coefficients K12 for OC
      call coagsub (r, imax, rcoag_ocn, rhob, rhooc, K12oc)   
      call coagsub (rp, imax, rcoag_ocn, rhob, rhooc, Kp12oc)   

c     Wavelength dependent complex rafractive indices (cref) 
c     are found from tabulated values for each aerosol component.
      call tabrefind (kcomp, ib, xlam, cref)

c     Aerosol hygroscopicities for max rh in the look-up tables (LUT),
c     for use as info in the header of each LUT:
      rh=relh(10)
      call hygro (rh, xbc, xdst, xoc, xs, xa, xss,
     $                      rhda, rhca, rhdss, rhcss)

c     Open output files for use in CAM-Oslo 
      call openfiles(kcomp,iopt,itilp,ib)

c     Adding header info for all look-up tables for each kcomp (only kcomp*.out yet):
csoa      call tableinfo ( kcomp, iSOA, xbc, xdst, xoc, xs, xa, xss,
csoa     $ relh, catote, catot, frac, frabc, fraq, ib, ibcam, itilp)
      call tableinfo ( kcomp, xbc, xdst, xoc, xs, xa, xss, relh, 
     $ frombg, frbcbg, catote, catot, frac, frabc, fraq, ib, ibcam, 
     $ itilp)


c          Editable input values    ! Full range for look-up tables
                isup1  = 1          ! 1 - 9    Note:
                isup2  = 9          ! 1 - 9    no loop if iopt=1

                irelh1 = 1          ! 1 - 10   Note:
                irelh2 = 10         ! 1 - 10   no loop if iopt=0

                ifombg1 = 1        ! 1 - 6    Note:
                ifombg2 = 6        ! 1 - 6    loop only for kcomp=1

                ifbcbg1 = 1        ! 1 - 6    Note:
                ifbcbg2 = 6        ! 1 - 6    loop only for kcomp=4

                ictot1 = 1          ! 1 - 6    Note: loop over       
                ictot2 = 6          ! 1 - 6    ictot OR ictote

                ictote1= 1          ! 1 - 16   Note: loop over ictot      
                ictote2= 16         ! 1 - 16   OR ictote, not both

                ifac1 = 1           ! 1 - 6    Note:
                ifac2 = 6           ! 1 - 6    no loop if kcomp=0

                ifbc1 = 1           ! 1 - 6    Note:
                ifbc2 = 6           ! 1 - 6    no loop if kcomp=0-4

                ifaq1 = 1           ! 1 - 6    Note:
                ifaq2 = 6           ! 1 - 6    no loop if kcomp=0-3

c          Do not edit the following input values!!! 
                if(iopt.eq.1) then   ! no supersaturation loop
                 isup1  = 1
                 isup2  = 1
                else                 ! no RH loop
                 irelh1 = 1
                 irelh2 = 1
                 if(itilp.eq.1) isup2 = 1 ! lognormal size fitting only for dry aerosols
                endif
                if(kcomp.ne.1.or.itilp.eq.1) then  !no OM (as SOA) internally mixed in the background
                 ifombg1 = 1
                 ifombg2 = 1
                endif
                if(kcomp.ne.4) then  !no BC internally mixed in the background
                 ifbcbg1 = 1
                 ifbcbg2 = 1
                endif
                if(kcomp.eq.0) then
                 ifac1  = 1
                 ifac2  = 1
                 ifbc1  = 1
                 ifbc2  = 1
                 ifaq1  = 1
                 ifaq2  = 1
                 ictote1= 1
                 ictote2= 1
                 ictot1 = 1   
                 ictot2 = 1
                elseif(kcomp.ge.1.and.kcomp.le.3) then
                 ifbc1  = 1
                 ifbc2  = 1
                 ifaq1  = 1
                 ifaq2  = 1
                 ictot1 = 1   
                 ictot2 = 1
csoa                elseif(kcomp.eq.4) then  ! background is OC, so that all added OC or BC comes
csoa                 ifbc1  = 1              ! as BC (fac here means BC/(BC+OC) for the background),
csoa                 ifbc2  = 1              ! and is homogeneously mixed (wrt. r).
csoa                 ictot1 = 1              ! Added SO4 is distributed     
csoa                 ictot2 = 1              ! according to D'(r) or r>rc, however. 
                elseif(kcomp.eq.4) then  ! background is OC and BC, and all added carbonaceous 
                 ifbc1  = 1              ! comes as SOA (fac=SOA/(SOA+Sulfate) added).
                 ifbc2  = 1              ! BC and OC is homogeneously mixed (wrt. r) in the 
                 ictot1 = 1              ! background. Added SO4 and SOA are distributed     
                 ictot2 = 1              ! according to D'(r) or r>rc (for sulfate), however. 
                else
                 ictote1= 1   
                 ictote2= 1 
                endif  ! kcomp
                if(itilp.eq.1) then  !no fombg or fbcbg dependency for itilp=1
                 ifombg1 = 1
                 ifombg2 = 1
                 ifbcbg1 = 1
                 ifbcbg2 = 1
                endif

                do 540 isup = isup1, isup2
c                do 540 isup = 1,1

              do 540 irelh = irelh1, irelh2
c               do 540 irelh = 1,10,9

c              supersaturation or relative humidity: the value of iopt determines which 
c              one (supers or rh) is used in the calculations.
               supers=sup(isup)
               rh=relh(irelh)
               if(iopt.eq.0) then
                 rh=0.05
                endif
cX              extra test loop for hygroscopic growth plots (with e.g. irelh=1,1 in the loop above)
c                do 540 irh=1,99,2
c                rh=0.01*real(irh)
c                or 
c                do 540 irh=1,199
c                rh=0.005*real(irh)
cX

c               Aerosol hygroscopicities (RH dependent) 
c               and points of deliquescence & crystallisation 
                call hygro (rh, xbc, xdst, xoc, xs, xa, xss,
     $                      rhda, rhca, rhdss, rhcss)

              do 540 ifombg = ifombg1, ifombg2
c              do 540 ifombg = 4,4

              do 540 ifbcbg = ifbcbg1, ifbcbg2
c              do 540 ifbcbg = 1,1

            do 540 ictot  = ictot1, ictot2
            do 540 ictote = ictote1, ictote2
c            do 540 ictot = 1,5,4
c            do 540 ictote = 1,10,9

          do 540 ifac = ifac1, ifac2
c          do 540 ifac = 1,5,4

         do 540 ifbc = ifbc1, ifbc2
c        do 540 ifbc = 1,5,4

      do 540 ifaq = ifaq1, ifaq2
c      do 540 ifaq = 1,5,4


      if(kcomp.eq.1) then
       write(*,*) 'kcomp,irelh,ifombg,ictote,ifac=',
     $ kcomp,irelh,ifombg,ictote,ifac
      elseif(kcomp.eq.1.or.kcomp.eq.3) then
       write(*,*) 'kcomp,irelh,ictote,ifac=',kcomp,irelh,ictote,ifac
      elseif(kcomp.eq.4) then
       write(*,*) 'kcomp,irelh,ifbcbg,ictote,ifac,ifaq=',
     $ kcomp,irelh,ifbcbg,ictote,ifac,ifaq
      else 
       write(*,*) 'kcomp,irelh,ictot,ifac,ifbc,ifaq=',
     $ kcomp,irelh,ictot,ifac,ifbc,ifaq 
      endif

c     Basic input parameters to the table calculations:

c     supersaturation or relative humidity: the value of iopt determines which 
c     one (supers or rh) is used in the calculations.
c      supers=sup(isup)
c      rh=relh(irelh)
c      if(iopt.eq.0) then
c        rh=0.05
c      endif

cX     extra test loop for hygroscopic growth (with e.g., irelh=1,1)
c      do 540 irh=1,99
c      rh=0.01*real(irh)
c      do 540 irh=1,199
c      rh=0.005*real(irh)
cX

c     concentrations of internally mixed SO4, BC and OC. Cas1, Cas2, Cas3 
cM     and Caso4 is internally mixed SO4 from condensation, coagulation, 
cM     cloud processing and all of the above, respectively. Cabc and Caoc is 
c     and Caso4 is internally mixed SO4 from condensation (H2SO4), coagulation
c     (H2SO4), cloud processing ((NH4)2SO4) and all of the above, respectively. 
c     Cabc and Caoc is all internally mixed BC and OC (from coagulation), respectively.
cSOA  Similarly, casoa is internally mixed OC from condensation.
csoa      if(itot.eq.0.or.(kcomp.le.4.and.ictote.eq.1)) then
csoa    initializing variables for concentrations of added mass onto the background
        Cas1=1.e-40 
        Cas2=1.e-40 
        Cas3=1.e-40
        Cabc=1.e-40
        Caoc=1.e-40
csoa      else
csoa        if(kcomp.ge.1.and.kcomp.le.3) then
        if(kcomp.ge.1.and.kcomp.le.4) then
          Cac=frac(ifac)*catote(ictote)         !  added Carbonaceous from condensation (SOA)
csoa        elseif(kcomp.eq.4) then
csoa          Cac=1.e-40
        else
          Cac=frac(ifac)*catot(ictot)           !  added Carbonaceous from condensation (SOA) and coagulation 
        endif
        Cabc=frabc(ifbc)*Cac                    !  added BC from coagulation 
        if(Cabc.lt.1.e-40) Cabc=1.e-40  
        Caoc=(1.0-frabc(ifbc))*Cac              !  added OC from condensation and coagulation
        if(Caoc.lt.1.e-40) Caoc=1.e-40
csoa        if(kcomp.ge.1.and.kcomp.le.3) then
        if(kcomp.ge.1.and.kcomp.le.4) then
          Caso4=(1.0-frac(ifac))*catote(ictote) !  added Sulfate from condensation (H2SO4)
csoa        elseif(kcomp.eq.4) then
csoa          Caso4=catote(ictote)                  !  added Sulfate from condensation (H2SO4) and wet phase ((NH4)2SO4)
        else
          Caso4=(1.0-frac(ifac))*catot(ictot)     !  added Sulfate from condensation and coagulation (H2SO4) and wet phase ((NH4)2SO4)
        endif
        if(Caso4.lt.1.e-40) Caso4=1.e-40        
        if(kcomp.ge.1.and.kcomp.le.4) then
          Cas1=(1.0-fraq(ifaq))*Caso4           !  added Sulfate from condensation (H2SO4)
          Cas2=1.e-40                           !  lump coagulation with condensation for these modes
        else
          Cas1=1.e-40                           !  lump condensation with coagulation for these modes
          Cas2=(1.0-fraq(ifaq))*Caso4           !  added Sulfate from coagulation (H2SO4)
        endif
        Cas3=fraq(ifaq)*Caso4                   !  added Sulfate from wet phase production ((NH4)2SO4)
        if(Cas1.lt.1.e-40) Cas1=1.e-40
        if(Cas2.lt.1.e-40) Cas2=1.e-40
        if(Cas3.lt.1.e-40) Cas3=1.e-40
csoa      endif

      Caso4=Cas1+Cas2+Cas3
cM      fr=Cas3/Caso4               ! wet-phase fraction of SO4
      faq=fraq(ifaq)                            !  wet-phase mass fraction of added sulfate (H2SO4 or (NH4)2SO4)
      fac=frac(ifac)                            !  Carbonaceous mass fraction of total added mass
csoa      fabc=frabc(ifbc)                          !  BC mass fraction of added Carbonaceous mass (or background BC/(BC+OC) for kcomp=4)
      fabc=frabc(ifbc)                          !  BC mass fraction of added Carbonaceous mass
csoa
      fombg=frombg(ifombg)
c     fombg is the OM (as SOA) mass fraction in the background SO4&SOA(Ait) mode.  
c     The respective volume fraction of OM in background is then:
      vombg=1.0/(1.0+(1.0-fombg)/(fombg*rhosv/rhooc+eps))
      fbcbg=frbcbg(ifbcbg)
c     fbcbg is the BC mass fraction in the background OC&BC(Ait) mode.  
c     The respective volume fraction of BC in background is then:
      vbcbg=1.0/(1.0+(1.0-fbcbg)/(fbcbg*rhooc/rhobc+eps))
csoa      

      if(kcomp.ge.1.and.kcomp.le.10) then
c       contribution to Ctot from the background mode
csoa       if(kcomp.eq.4) then
       if(kcomp.eq.1) then
        Ctot0=Ctotnull*(1.0+vombg*(rhooc/rhob-1.0)) ! -> Ctotnull*0.815 for ren OM (vombg=fombg=1).
       elseif(kcomp.eq.4) then
csoa        Ctot0=Ctotnull*(fac*rhobc/rhooc+1.0-fac)    ! -> Ctotnull*1.333 for ren BC (fac=1).
        Ctot0=Ctotnull*(1.0+vbcbg*(rhobc/rhob-1.0)) ! -> Ctotnull*1.333 for ren BC (vbcbg=fbcbg=1).
        write(*,*) 'Ctotnull =', Ctot0
       endif
        write(999,*) 'background contribution:'
        write(999,*) Ctot0
cM       contribution to Ctot from internally mixed SO4
c       contribution to Ctot from internally mixed (non-background) H2SO4 and (NH4)2SO4
c       (note: only H2SO4 for kcomp=1 since ifaq=1 there)
cM        dCtot=(fr*(Ms/Mso4)+(1.0-fr)*(Msv/Mso4))*Caso4  
        dCtot=Caso4 
        Ctot=Ctot0+dCtot
        write(999,*) 'sulfate contribution (a, tot):'
        write(999,*) dCtot, Ctot
c       contribution to Ctot from internally mixed (non-background) BC
        dCtot=Cabc
        Ctot=Ctot+dCtot
        write(999,*) 'bc contribution (a, tot):'
        write(999,*) dCtot, Ctot
c       contribution to Ctot from internally mixed (non-background) OC
        dCtot=Caoc
        Ctot=Ctot+dCtot
        write(999,*) 'oc contribution (a, tot):'
        write(999,*) dCtot, Ctot
      else
        Ctot=Ctot0
      endif
      write(*,*) 'dry Ctot =', Ctot

      if(kcomp.ge.1.and.kcomp.le.4) then
        cat=catote(ictote)
      else
        cat=catot(ictot)
      endif

c     Calculate modified dry size distributions for process specific 
c     SO4 and BC (and OC) internally mixed with the background aerosol
ccccc6ccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      call conteq (r, rp, rbcn, d, itot, imax, ictot, ictote, ifaq, 
cM     $ imini, Ms, Msv, Mso4, rhos, rhobc, rhooc, rhob, rhosv, rhoc2,
     $ imini, rhos, rhobc, rhooc, rhob, rhosv, rhoc2,
     $ Nnatk, fcondk, fcoagk, faqk, Cas1, Cas2, Cas3, Cabc, Caoc, Ctot0,
     $ dndlrk0, dndlrkny, ntot, Dmsoa, Dmpsoa, Dm, Dmp, K12, Kp12, 
     $ K12oc, Kp12oc, K12so4, Kp12so4, ismolar, vbci, voci, vsi, vai, 
     $ cintbg, cintsc, cintsa, cintbc, cintoc, cintbg05, cintsc05, 
     $ cintsa05, cintbc05, cintoc05, cintbg125, cintsc125, cintsa125, 
     $ cintbc125, cintoc125, aaero, aaeros, vaero, vaeros, fracdim, 
csoa     $ kcomp, fac, iSOA)
     $ kcomp, vombg, vbcbg, fac)
cSOA

c     Hygroscopic growth is taken into account in subroutine rhsub,
c     either for the given relative humidity (if iopt=1), or for the
c     given supersaturation (if iopt=0). CCN results are written to
c     file from rhsub. 
      Cdry=Ctot
      if(ksol.eq.1.and.itilp.eq.0) then
        call rhsub (imax, rh, d, r, rp, dndlrkny, vsi, vbci, voci, 
csoa     $   vssol, vbcsol, vocsol, vasol, vw, fki, itot, rhos, 
     $   fombg, fbcbg, vombg, vbcbg, 
     $   vssol, vbcsol, vocsol, vasol, vw, fki, itot, rhos, 
     $   rhosv, rhobc, rhooc, rhob, rhow, Ctot, kcomp, iopt, supers, 
csoa     $   rsup, ismolarh, iSOA, cat, fac, fabc, faq, CCN,
     $   rsup, ismolarh, cat, fac, fabc, faq, CCN,
     $   xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)
      endif

c     Tabulate aerosol size distribution after hygroscopic growth, 
c     and check how well total aerosol number is conserved
      numb=0.0
      if(itot.eq.0) then
        do i=1,imax
          dndlrk(i)=dndlrkny(i)
          numb=numb+dndlrk(i)*d
ct
c          write(12,100) r(i), dndlrk(i) 
c          write(14,100) r(i), dndlrk(i)*(4.0*pi/3.0)*r(i)**3
ct
          if(iopt.eq.1)
     $     write(9001,400) r(i), dndlrk(i), rh, kcomp 
        enddo
      else
        do i=1,imax
          if(dndlrkny(i).lt.0.0) then
            write(*,*) 'dndlrkny(i) < 0 !'
            stop
          endif
          numb=numb+dndlrkny(i)*d
c          write(13,100) r(i), dndlrkny(i)
c          write(14,100) r(i), dndlrkny(i)*(4.0*pi/3.0)*r(i)**3
c          if(iopt.eq.1.or.itilp.eq.1) then
          if((iopt.eq.1.or.itilp.eq.1).and.ib.ne.19) then
           write(9001,500) r(i), dndlrkny(i),  
     $      cat, fac, fabc, faq, rh, kcomp
          endif
        enddo
        if(itilp.eq.1) then
csoa     Note tha dry lognormal the fitted size parameters do not depend
csoa     on the mass fraction fombc for kcomp=1     ! or fac for kcomp=4 SJEKK!!!!!!!!!           
         call modetilp(pi, imax, d, r, dndlrkny, dndlrk0,
csoa     $     cat, fac, fabc, faq, kcomp, iSOA)
     $     cat, fac, fabc, faq, kcomp)
        endif
      endif
c      write(*,*) 'numbny=', numb

c     Sizemie determines the spectral aerosol gross (size integrated) 
c     optical parameters (by calling the Mie code for each particle size), 
c     and writes the result to file.
      if(iopt.eq.1) then 
csoa      call sizemie(imini, imaxi, iSOA, r, rbcn, d, vsi, vbci, voci, vai, 
      call sizemie(imini, imaxi, r, rbcn, d, vsi, vbci, voci, vai,
     $  vombg, fombg, vbcbg, fbcbg,                                          
     $  dndlrk, dndlrkny, kcomp, itot, ib, vssol, vbcsol, vocsol, vasol, 
     $  vw, fki, rh, Ctot, Nnat, cat, fac, fabc, faq, fracdim, xlam, 
     $  xlami, xlamb, xlame, fband, fb, cref, omega, gass, bext, kext)
      endif

      if(iopt.eq.1) then
c        write(*,*) 'rh, Caso4, Cas1, Cas3, Cabc'
c        write(*,300) rh, Caso4, Cas1, Cas3, Cabc
c        write(*,*)
c        write(*,*) 'irelh, ictot, ifbc, ifaq'
c        write(*,*) irelh, ictot, ifbc, ifaq
c        write(1160+kcomp,*) rh, omega(9)
c        write(1170+kcomp,*) rh, gass(9)
c        write(1180+kcomp,*) rh, bext(9)
c        write(1190+kcomp,*) rh, kext(9)
c
c      Here comes the aerodryk*.out look-up tables:
c
      if(ib.eq.31.and.irelh.eq.1) then
       aaerol=aaero-aaeros                              
       vaerol=vaero-vaeros
       if(itot.eq.1) then                              
         if(kcomp.eq.1) then
csoa           write(9600,2000) kcomp, cat, fac,
           write(9600,2100) kcomp, fombg, cat, fac,
     $     cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  
     $     cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,        
     $     cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol          
csoa       write(*,*) 'cintoc =', cintoc
         elseif(kcomp.eq.2.or.kcomp.eq.3) then         
csoa           write(9600,2000) kcomp, cat,
           write(9600,2000) kcomp, cat, fac,
     $     cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  
     $     cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,        
     $     cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol          
         elseif(kcomp.eq.4) then
csoa           write(9600,2100) kcomp, cat, fac, faq,
           write(9600,3000) kcomp, fbcbg, cat, fac, faq,
     $     cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  
     $     cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,        
     $     cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol          
         else  ! (kcomp=5-10))
           write(9600,3000) kcomp, cat, fac, fabc, faq, 
     $     cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  
     $     cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,        
     $     cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol          
         endif 
       else  ! itot (kcomp=0)
         write(9600,4000) kcomp, cintbg, cintbg05, cintbg125, 
     $     aaeros, aaerol, vaeros, vaerol          
       endif ! itot
      endif ! ib & relh
c
      endif  ! iopt

 540  continue         ! ifaq, ifbc, ifac, ictot/ictote, ifombg, ifbcbg, irelh/isup

      close(9000) 
      close(9001)
      close(9002) 
      close(9003) 
      close(9600) 

      enddo  ! kcomp
     

 100  format(2(x,e10.4))
csoa 200  format(A25,6I3)
 300  format(x,f9.3,5(x,e9.3))
 400  format(2(x,e12.5),f7.2,I3)
 500  format(6(x,e12.5),f7.2,I3)
csoa 2000 format(I2,e10.3,19e10.3)
csoa 2100 format(I2,2e10.3,19e10.3)
 2000 format(I2,21e10.3)
 2100 format(I2,22e10.3)
csoa 2500 format(I2,3e10.3,19e10.3)
csoa 3000 format(I2,3e10.3,f5.2,19e10.3)
 3000 format(I2,23e10.3)
 4000 format(I2,7e11.4)
       

      end
