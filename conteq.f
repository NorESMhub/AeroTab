      subroutine conteq (r, rp, rbcn, d, itot, imax, ictot, 
cM     $ ictote, ifaq, imini, Ms, Msv, Mso4, rhos, rhobc, rhooc, rhob, 
     $ ictote, ifaq, imini, rhos, rhobc, rhooc, rhob, 
     $ rhosv, rhoc2, Nnatk, fcondk, fcoagk, faqk, Cas1, Cas2, Cas3, 
     $ Cabc, Caoc, Ctot0, dndlrk, dndlrkny, ntot, Dmsoa, Dmpsoa, Dm,
     $ Dmp, K12in, Kp12in, K12ocin, Kp12ocin, K12so4in, Kp12so4in, 
     $ ismolar, vbci, voci, vsi, vai, cintbg, cintsc, cintsa, cintbc, 
     $ cintoc, cintbg05, cintsc05, cintsa05, cintbc05, cintoc05, 
     $ cintbg125, cintsc125, cintsa125, cintbc125, cintoc125, aaero, 
csoa     $ aaeros, vaero, vaeros, fracdim, kcomp, fac, iSOA) 
     $ aaeros, vaero, vaeros, fracdim, kcomp, vombg, vbcbg, fac) 
cSOA

c **********************************************************************************
c     Created by Alf Kirkevåg.
c **********************************************************************************

c     Here the modified dry size distributions for process specific 
c     SO4, BC and OC internally mixed with the background aerosol is 
c     calculated. 
c     NOTE: For kcomp=4 (the OC&BC(Ait) mode) we assume that the
c     background mode consists of OC, and add BC homogeneously wrt.
c     radius (except at very small radii due to small OC amounts due
c     to a frac a little <1), if no so4 from condensation is added.
c     Therefore we assume that dvbc=0 in calculating modified size
c     distributions, but take the homogeneously internally mixed BC
c     into account afterwards, by using an adjusted cbc(i) wrt. the
c     volume fractions vbci etc. (for use in rhsub and sizemie) as 
c     well as by remembering to multiply the normalized N4 numbers 
c     by a constant when the tables are to be used in CAM-Oslo.  
csoa
c     New treatment for kcomp=4: fac is now mass fraction of added OC
c     (as SOA), as for the other modes, while instead fbcbg is the mass 
c     fraction of BC in the background OC&BC(Ait) mode. Similarly, fombg 
c     is the mass fraction of OM in the background SO4&SOA(ait) mode. 
csoa 
ccccc6ccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      implicit none

      INTEGER  i, ic, imax, ix, imini, itot, ictot, ictote, 
     $ ismolar, ifaq, j, jmax, jmaxx, k, kcomp, itest
      REAL K12in(0:101), Kp12in(0:101), K12ocin(0:101), 
     $ Kp12ocin(0:101),  K12so4in(0:101), Kp12so4in(0:101)
      REAL r(0:100), rp(0:100), dip(0:100), Dm(0:100), 
     $ Dmp(0:100), K12(0:101), Kp12(0:101), dninc(0:100), 
     $ dndlrk(0:100), dndlrkny(0:100), fracdim(0:100), rhorbc(0:100),
     $ vcbg(100), vcbc(100), vcoc(100), vcsu12(100), vcsu3(100), 
     $ vcsu(100), vbci(0:100), voci(0:100), vsi(0:100), vai(0:100), 
     $ dqsu12(100), dqsu3(100), dqbc(100), dqoc(100), cbg(0:100), 
     $ csu12(0:100), csu3(0:100), csu(0:100), cbc(0:100), coc(0:100), 
     $ dcincbg(100), dcincs12(100), dcincs3(100), dcincbc(100), 
     $ dcincoc(0:100), dncny(0:100), K12oc(0:101), Kp12oc(0:101),
     $ K12so4(0:101), Kp12so4(0:101)
      REAL rbcn, rc, rcmin, rcmax,  rjm, rjmg, rjmax, rjmaxx, d,
     $ Nnatk, ntot, nt, Nag, NrD, NK12, NK12oc, NK12so4, fcondk, fcoagk, 
     $ faqk, fr, frcoag, radikand, dv, dvcon, dvcos, dvcoa, dvcoaoc, 
     $ dvaq, dvbc, dvoc, dvaq0, dvs1, dvs2, dvs3, rhos, rhooc, rhobc, 
cM     $ rhob, rhosv, rhoc2, Ms, Msv, Mso4, Cas1, Cas2, 
     $ rhob, rhosv, rhoc2, Cas1, Cas2, 
     $ Cas3, Caso4, Cabc, Caoc, Ctot0, cintbg, cintsu, cintsc, cintsa, 
     $ cintbc, cintoc, cintbg05, cintsu05, cintsc05, cintsa05, 
     $ cintbc05, cintoc05, cintbg125, cintsc125, cintsa125, cintbc125, 
     $ cintoc125, dcintbg, vtot, aaero, aaeros, vaero, vaeros, 
     $ pi, e, p1, p2, fac, vaisave(0:100) 
cSOA  Added May & December 2013
      REAL Dmsoa(0:100), Dmpsoa(0:100), NrDsoa, dvsoa, dvconsoa 
csoa  July 2015
      REAL vombg, vbcbg, rhobg
csoa      INTEGER iSOA
cSOA

c     critical radius for cloud processing, rc, ranging between rcmin and 
c     rcmax, from Chuang and Penner (1995). We have assumed that this
c     range is independent of background aerosol 
      PARAMETER (rcmin=0.05, rcmax=0.2)

      PARAMETER (pi=3.141592654, e=2.718281828 )
csoa      PARAMETER (eps=1.e-50)

cM      Caso4=Cas1+Cas2+Cas3
cM      frcoag=Cas2/Caso4 
cM      fr=Cas3/Caso4
      Caso4=Cas1+Cas2+Cas3        ! total mass conc. of H2SO4 and (NH4)2SO4
      frcoag=Cas2/Caso4           ! (H2SO4 coagulate)/Caso4 
      fr=Cas3/Caso4               ! (wet-phase (NH4)2SO4)/Caso4
c      write(*,*) 'Cas1,2,3=', Cas1, Cas2, Cas3
c      write(*,*) 'fr, frcoag=', fr, frcoag

c     Initial guess of jmax, the maximum required iterations to satisfy
c     the stability criterion for the continuity equation. 
      jmaxx=10000 

c     Initiate modified size distribution calculations. Estimate
c     amount of "moves to the right" (ix) so that jmax < jmaxx,
c     and search for a sufficiently large total iteration number, 
c     jmax, to satisfy the stability criterium.

c     Initially modified size distribution = background size distribution 

      do i=0,imax
         K12(i)    = K12in(i)
         Kp12(i)   = Kp12in(i)
         K12oc(i)  = K12ocin(i)
         Kp12oc(i) = Kp12ocin(i)
         K12so4(i) = K12so4in(i)
         Kp12so4(i)= Kp12so4in(i)
c       write(19,*) r(i), Dm(i)
c       write(20,*) r(i), Dmsoa(i)
      enddo

      itest=0
 11   do i=0,imax                     
        dndlrkny(i)=dndlrk(i)
      enddo 
      if(itot.eq.0) then
        ix=0                  
      else
        if(itest.eq.0) ix=1  
      endif
 12   do i=0,imax
       if(i.le.ix) then
         K12(i)=0.0
         Kp12(i)=0.0
         K12oc(i)=0.0
         Kp12oc(i)=0.0
         K12so4(i)=0.0
         Kp12so4(i)=0.0
         dndlrkny(i)=1e-50
       endif
      enddo

c     Below, condensation of H2SO4 or SOA and coagulation of BC, OC and SO4 aerosol
c     onto the background distribution for the first time step is calculated   
cSOA  For SOA, add condensation --> NrDsoa (NrD is for so4 only)!
      NrD=0.0             ! H2SO4
      NrDsoa=0.0          ! SOA
      NK12=0.0            ! BC
      NK12oc=0.0          ! OC (OM)
      NK12so4=0.0         ! H2SO4
      do i=0,imax
        NrD=NrD+dndlrkny(i)*Dm(i)*r(i)*d
        NrDsoa=NrDsoa+dndlrkny(i)*Dmsoa(i)*r(i)*d  ! SOA
        NK12=NK12+dndlrkny(i)*K12(i)*d
        NK12oc=NK12oc+dndlrkny(i)*K12oc(i)*d
        NK12so4=NK12so4+dndlrkny(i)*K12so4(i)*d
      enddo

c     Process specific volumes per volume of dry air to be added:     
c     cloud processed sulfate (s3) is assumed to exist as (NH4)2SO4,
c     while sulfate from diffusional growth (s1) and coagulation (s2) 
c     exist (just like the nucleation mode) as H2SO4
cM      dvs1=(fcondk/Nnatk)*1e-9*(1.0-frcoag-fr)*Caso4*(Msv/Mso4)/rhosv  
cM      dvs2=(fcoagk/Nnatk)*1e-9*frcoag*Caso4*(Msv/Mso4)/rhosv         
      dvs1=(fcondk/Nnatk)*1e-9*(1.0-frcoag-fr)*Caso4/rhosv  ! volume of H2SO4 condensate
      dvs2=(fcoagk/Nnatk)*1e-9*frcoag*Caso4/rhosv           ! volume of H2SO4 coagulate
      dvbc=(fcoagk/Nnatk)*1e-9*Cabc/rhobc                   ! volume of BC coagulate 
cSOA      dvoc=(fcoagk/Nnatk)*1e-9*Caoc/rhooc                    
cSOA  assuming fcondsoa/Nnatk=1 as above (/=1 is never used) 
csoa      if(kcomp.eq.1.and.iSOA.eq.1) then
csoa      if(kcomp.ge.1.and.kcomp.le.3) then
      if(kcomp.ge.1.and.kcomp.le.4) then        ! OC only comes as SOA 
        dvsoa=1e-9*Caoc/rhooc
        dvoc=1.e-50
      else                                      ! SOA is lumped together with and treated as OC coagulate 
        dvsoa=1.e-50
        dvoc=(fcoagk/Nnatk)*1e-9*Caoc/rhooc
      endif
cSOA

c     searching for sufficiently large total number of iterations, 
c     jmax (>10000), to satisfy the stability criterium for the 
c     continuity equation
      rjmg=0.0
      do i=ix+1,imax
        dvcon=dvs1*r(i)*Dm(i)/NrD
        dvcos=dvs2*K12so4(i)/NK12so4
        dvcoa=dvbc*K12(i)/NK12
        dvcoaoc=dvoc*K12oc(i)/NK12oc
cSOA
csoa        if(kcomp.eq.1.and.iSOA.eq.1) then
        dvconsoa=dvsoa*r(i)*Dmsoa(i)/NrDsoa
csoa          dv=dvcon+dvcoa+dvcos+dvcoaoc+dvconsoa   
csoa        else
csoa          dvconsoa=0.0
csoa          dv=dvcon+dvcoa+dvcos+dvcoaoc   
csoa        endif
        dv=dvcon+dvcoa+dvcos+dvcoaoc+dvconsoa  ! csoa 
cSOA
        rjm=3e12*dv/(4.0*pi*r(i)**3.0*((1.0+d/log10(e))**3.0-1.0))
        rjm=rjm/(1.01-fr)
        if(rjm.gt.rjmg) then
          rjmax=rjm
c          write(*,*) i, r(i), rjmax 
        endif
        rjmg=rjm
      enddo

      rjmaxx=1.0*jmaxx
      if(rjmax.gt.rjmaxx) then
        ix=ix+1
        if(ix.gt.imini) then
ct          ix=1
          itest=1
          jmaxx=jmaxx+jmaxx
c          write(*,*) jmaxx
          goto 11
        endif
        goto 12
      endif

      jmax=int(rjmax)+1
      if(jmax.lt.10000.and.(ictot.gt.1.or.ictote.gt.1)) then
        jmax=10000
        if(ifaq.eq.6) jmax=20000
      endif
c      write(*,*) 'jmax, ix =', jmax, ix 
c      write(*,*) 'Caoc =', Caoc
     
c     Process specific volumes of SO4, BC and OC per volume of dry air 
c     to be added PER ITERATION is then determined:     
cM      dvs1=(fcondk/Nnatk)*1e-9*(1.0-frcoag-fr)*Caso4*(Msv/Mso4)
cM     $     /(rhosv*jmax)                                        
cM      dvs2=(fcoagk/Nnatk)*1e-9*frcoag*Caso4*(Msv/Mso4)/(rhosv*jmax)
cM      dvs3=(faqk/Nnatk)*1e-9*fr*Caso4*(Ms/Mso4)/(rhos*jmax)         
      dvs1=(fcondk/Nnatk)*1e-9*(1.0-frcoag-fr)*Caso4
     $     /(rhosv*jmax)                                        
      dvs2=(fcoagk/Nnatk)*1e-9*frcoag*Caso4/(rhosv*jmax)
      dvs3=(faqk/Nnatk)*1e-9*fr*Caso4/(rhos*jmax)         
      dvbc=(fcoagk/Nnatk)*1e-9*Cabc/(rhobc*jmax)                           
cSOA      dvoc=(fcoagk/Nnatk)*1e-9*Caoc/(rhooc*jmax)                           
csoa      if(kcomp.eq.1.and.iSOA.eq.1) then
      if(kcomp.ge.1.and.kcomp.le.4) then
        dvsoa=1e-9*Caoc/(rhooc*jmax)
        dvoc=1.e-50
      else
        dvsoa=1.e-50
        dvoc=(fcoagk/Nnatk)*1e-9*Caoc/(rhooc*jmax)
      endif
cSOA

c     initialize arrays for mass concentrations of the 
c     background aerosol (after correcting for internal 
c     mixtures in the background aerosol)
csoa
      if(kcomp.eq.1) then
        rhobg=rhob*(1.0+vombg*(rhooc/rhob-1.0))
      elseif(kcomp.eq.4) then
        rhobg=rhob*(1.0+vbcbg*(rhobc/rhob-1.0))
      else
        rhobg=rhob
      endif
      do i=1,imax
csoa        cbg(i)=1.0e-3*(4.0*pi/3.0)*r(i)**3.0*(rhob*dndlrkny(i))  ! ug/m^3
        cbg(i)=1.0e-3*(4.0*pi/3.0)*r(i)**3.0*(rhobg*dndlrkny(i))  ! ug/m^3
      enddo
csoa
c     ... and initialize arrays for internally mixed 
cM     BC, sulfate from condensation and coagulation, sulfate from 
cM     cloud processing, and for the background aerosol   
csoa     BC, H2SO4 from condensation and coagulation, (NH4)2SO4 from 
c     BC, H2SO4 and OC (POM or SOA) from condensation and coagulation 
c     and (NH4)2SO4 from cloud processing
      do i=1,imax
        cbc(i)=1.0e-100
        coc(i)=1.0e-100
        csu12(i)=1.0e-100
        csu3(i)=1.0e-100                          
      enddo
csoa
c     take into account that the background is an internal mixture
c      if(kcomp.eq.1) then
c        do i=1,imax 
c          cbg(i)=cbg(i)*(1.0+vombg*(rhooc/rhob-1.0))             ! ug/m^3
c        enddo 
c        write(*,*) vombg
c      elseif(kcomp.eq.4) then
c        do i=1,imax 
c          cbg(i)=cbg(i)*(1.0+vbcbg*(rhobc/rhob-1.0))             ! ug/m^3
c        enddo 
c        write(*,*) vbcbg
c      endif
csoa

c     then solve continuity equation using jmax time steps/iterations 
      do 20 j=1,jmax

        rc=rcmin+j*(rcmax-rcmin)/jmax

c       initialization of key variables for each time step
        NrD=0.0
        NrDsoa=0.0          ! SOA
        NK12=0.0
        NK12oc=0.0
        NK12so4=0.0
        Nag=0.0
        k=0
c       variables for growth by condensation and coagulation
        do i=1,imax
          if(i.le.ix) dndlrkny(i)=1.0e-50 
          NrD=NrD+dndlrkny(i)*Dmp(i)*rp(i)*d
          NrDsoa=NrDsoa+dndlrkny(i)*Dmpsoa(i)*rp(i)*d  ! SOA
          NK12=NK12+dndlrkny(i)*Kp12(i)*d
          NK12oc=NK12oc+dndlrkny(i)*Kp12oc(i)*d
          NK12so4=NK12so4+dndlrkny(i)*Kp12so4(i)*d
          if(rp(i).ge.rc) k=k+1        
          if(k.eq.1) ic=i
        enddo
c       variables for growth by cloud processing (wetphase chemistry) 
        Nag=dndlrkny(ic)*log10(rp(ic)/rc)
        do i=ic+1,imax
          Nag=Nag+dndlrkny(i)*d
        enddo
        dvaq0=dvs3/Nag      ! as (NH4)2SO4

c       calculate process specific volumes of SO4 aerosol, BC and OC 
c       per volume of dry air to be added (per particle) in each size bin
        do i=1,imax       
          if(i.lt.ic) then
            dvaq=0.0
          elseif(i.eq.ic) then
            dvaq=dvaq0*(log10(rp(ic)/rc))/d
          elseif(i.ge.ic+1) then
            dvaq=dvaq0
          endif
          dvcon=dvs1*rp(i)*Dmp(i)/NrD        ! as H2SO4 
          dvcos=dvs2*Kp12so4(i)/NK12so4      ! as H2SO4
          dvcoa=dvbc*Kp12(i)/NK12  
          dvcoaoc=dvoc*Kp12oc(i)/NK12oc      ! NB
cSOA
csoa          if(kcomp.eq.1.and.iSOA.eq.1) then
          dvconsoa=dvsoa*r(i)*Dmsoa(i)/NrDsoa
csoa            dv=dvcon+dvaq+dvcoa+dvcos+dvcoaoc+dvconsoa   
csoa          else
csoa            dvconsoa=0.0
csoa            dv=dvcon+dvaq+dvcoa+dvcos+dvcoaoc  
csoa          endif
          dv=dvcon+dvaq+dvcoa+dvcos+dvcoaoc+dvconsoa  ! csoa 
cSOA
c         find the incement of log(r/um) at r=rp, i.e. 
c         in the center of the size bin, dip
          radikand=1.0+3.0e12*dv/(4.0*pi*rp(i)**3.0)
          dip(i)=log10(e)*(radikand**(1/3.0)-1.0)
c         process specific mass concentration increments (ug/m^3)
          dqbc(i)=1e9*rhobc*dvcoa*dndlrkny(i)
cSOA
csoa          if(kcomp.eq.1.and.iSOA.eq.1) then
          dqoc(i)=1e9*rhooc*(dvcoaoc+dvconsoa)*dndlrkny(i)
csoa          else
csoa            dqoc(i)=1e9*rhooc*dvcoaoc*dndlrkny(i)
csoa          endif
cSOA
          dqsu12(i)=1e9*rhosv*(dvcon+dvcos)*dndlrkny(i)  ! as H2SO4
          dqsu3(i)=1e9*rhos*dvaq*dndlrkny(i)             ! as (NH4)2SO4
        enddo

c       solve the continuity equations (using a simple upwind advection 
c       scheme) for the size distribution, dndlrkny, and for the process 
c       specific mass concentrations, dcinc*
        dip(0)=0.0
        do i=1,imax       
          dninc(i)=-(dndlrkny(i)*dip(i)-dndlrkny(i-1)*dip(i-1))/d
          dcincbg(i)=-(cbg(i)*dip(i)-cbg(i-1)*dip(i-1))/d
         if(ismolar.eq.0) then
          dcincbc(i)=-(cbc(i)*dip(i)-cbc(i-1)*dip(i-1))/d+dqbc(i)
          dcincoc(i)=-(coc(i)*dip(i)-coc(i-1)*dip(i-1))/d+dqoc(i)
          dcincs12(i)=-(csu12(i)*dip(i)-csu12(i-1)*dip(i-1))/d
     $                +dqsu12(i)
          dcincs3(i)=-(csu3(i)*dip(i)-csu3(i-1)*dip(i-1))/d+dqsu3(i)
         else
          dcincbc(i)=-(cbc(i)*dip(i)-cbc(i-1)*dip(i-1))/d
          dcincoc(i)=-(coc(i)*dip(i)-coc(i-1)*dip(i-1))/d
          dcincs12(i)=-(csu12(i)*dip(i)-csu12(i-1)*dip(i-1))/d
          dcincs3(i)=-(csu3(i)*dip(i)-csu3(i-1)*dip(i-1))/d
         endif
        enddo
        do i=1,imax       
          dndlrkny(i)=dndlrkny(i)+dninc(i)
          if(dndlrkny(i).lt.1.e-99) dndlrkny(i)=1.e-99 
          cbg(i)=cbg(i)+dcincbg(i)          
          coc(i)=coc(i)+dcincoc(i)          
          cbc(i)=cbc(i)+dcincbc(i)          
          csu12(i)=csu12(i)+dcincs12(i)          
          csu3(i)=csu3(i)+dcincs3(i)          
          csu(i)=csu12(i)+csu3(i)          
c          write(15,*) r(i), dndlrkny(i)
        enddo

c       here the anti-diffusive part of the upwind scheme by 
c       Smolarkiewicz (1983) kicks in, providing that the  
c       number of corrective steps is chosen larger than 0 
        if(ismolar.gt.0) then
c         size distribution (number concentration)
          call smolar (ismolar, imax, d, dndlrkny, dip)  
c          do i=1,imax       
c            write(16,*) r(i), dndlrkny(i)
c          enddo
c         background mass concentration
          call smolar (ismolar, imax, d, cbg, dip)  
c         non-backgrond OC mass concentration
          do i=1,imax
            dncny(i)=coc(i)
          enddo
          call smolar (ismolar, imax, d, dncny, dip)  
          do i=1,imax
            coc(i)=dncny(i)+dqoc(i)
          enddo
c         non-backgrond BC mass concentration
          do i=1,imax
            dncny(i)=cbc(i)
          enddo
          call smolar (ismolar, imax, d, dncny, dip)  
          do i=1,imax
            cbc(i)=dncny(i)+dqbc(i)
          enddo
cM         process specific and total SO4 mass concentrations
c         process specific and total (non-backgrond) H2SO4 or/and (NH4)2SO4 mass concentrations
          do i=1,imax
            dncny(i)=csu12(i)    ! as H2SO4
          enddo
          call smolar (ismolar, imax, d, dncny, dip)  
          do i=1,imax
            csu12(i)=dncny(i)+dqsu12(i)
          enddo
          do i=1,imax
            dncny(i)=csu3(i)     ! as (NH4)2SO4
          enddo
          call smolar (ismolar, imax, d, dncny, dip)  
          do i=1,imax
            csu3(i)=dncny(i)+dqsu3(i) 
            csu(i)=csu12(i)+csu3(i)   ! as H2SO4 + (NH4)2SO4 mass       
          enddo
        endif  ! ismolar

 20   continue    ! j=1,jmax 

c     check if total dry aerosol number is conserved
c     and calculate aerosol area and volume, total and below 0.5um
c     (thereby implicitely also above 0.5um), for AEROCOM diagnostics
      nt=0.0
      aaero=0.0
      vaero=0.0
      p1=4.0*pi
      p2=p1/3.0
      aaeros=0.99*p1*r(28)**2*dndlrkny(28)*d
      vaeros=0.99*p2*r(28)**3*dndlrkny(28)*d
      do i=1,imax       
        nt=nt+dndlrkny(i)*d
        aaero=aaero+p1*r(i)**2*dndlrkny(i)*d
        vaero=vaero+p2*r(i)**3*dndlrkny(i)*d
        if(i.le.27) then
          aaeros=aaeros+p1*r(i)**2*dndlrkny(i)*d
          vaeros=vaeros+p2*r(i)**3*dndlrkny(i)*d
        endif
      enddo  
c      write(*,*) 'Ntot og Nt er:', ntot, nt

c     size-integrated dry mass concentrations, integrated over all r,
c     and r<0.5um and r>1.25um (for AEROCOM).         
      cintbg=0.0
      cintbc=0.0
      cintoc=0.0
c      cintsu=0.0
      cintsc=0.0
      cintsa=0.0
      cintbg05=0.99*cbg(28)*d
      cintbc05=0.99*cbc(28)*d
      cintoc05=0.99*coc(28)*d
c      cintsu05=0.99*Mso4*(csu12(28)/Msv+csu3(28)/Ms)*d
cM      cintsc05=0.99*Mso4*(csu12(28)/Msv)*d
cM      cintsa05=0.99*Mso4*(csu3(28)/Ms)*d
      cintsc05=0.99*csu12(28)*d  ! as H2SO4
      cintsa05=0.99*csu3(28)*d   ! as (NH4)2SO4
      cintbg125=0.03*cbg(31)*d
      cintbc125=0.03*cbc(31)*d
      cintoc125=0.03*coc(31)*d
c      cintsu125=0.03*Mso4*(csu12(31)/Msv+csu3(28)/Ms)*d
cM      cintsc125=0.03*Mso4*(csu12(31)/Msv)*d
cM      cintsa125=0.03*Mso4*(csu3(31)/Ms)*d
      cintsc125=0.03*csu12(31)*d  ! as H2SO4
      cintsa125=0.03*csu3(31)*d   ! as (NH4)2SO4
      do i=1,imax
        if(cbg(i).lt.1.e-100)   cbg(i)=1.e-100
        if(cbc(i).lt.1.e-100)   cbc(i)=1.e-100
        if(coc(i).lt.1.e-100)   coc(i)=1.e-100
        if(csu12(i).lt.1.e-100) csu12(i)=1.e-100
        if(csu3(i).lt.1.e-100)  csu3(i)=1.e-100
        csu(i)=csu12(i)+csu3(i)  ! as H2SO4 + (NH4)2SO4 mass        
        cintbg=cintbg+cbg(i)*d
        cintbc=cintbc+cbc(i)*d
        cintoc=cintoc+coc(i)*d
c        cintsu=cintsu+(csu12(i)+csu3(i))*d
        cintsc=cintsc+csu12(i)*d  ! as H2SO4
        cintsa=cintsa+csu3(i)*d   ! as (NH4)2SO4
        if(i.le.27) then
          cintbg05=cintbg05+cbg(i)*d
          cintbc05=cintbc05+cbc(i)*d
          cintoc05=cintoc05+coc(i)*d
c          cintsu05=cintsu05+Mso4*(csu12(i)/Msv+csu3(i)/Ms)*d
cM          cintsc05=cintsc05+csu12(i)*d*Mso4/Msv
cM          cintsa05=cintsa05+csu3(i)*d*Mso4/Ms
          cintsc05=cintsc05+csu12(i)*d  ! as H2SO4
          cintsa05=cintsa05+csu3(i)*d   ! as (NH4)2SO4
        endif
        if(i.ge.32) then
          cintbg125=cintbg125+cbg(i)*d
          cintbc125=cintbc125+cbc(i)*d
          cintoc125=cintoc125+coc(i)*d
c          cintsu125=cintsu125+Mso4*(csu12(i)/Msv+csu3(i)/Ms)*d
cM          cintsc125=cintsc125+csu12(i)*d*Mso4/Msv
cM          cintsa125=cintsa125+csu3(i)*d*Mso4/Ms
          cintsc125=cintsc125+csu12(i)*d  ! as H2SO4
          cintsa125=cintsa125+csu3(i)*d   ! as (NH4)2SO4
        endif
      enddo
c****************************************
cM      if(kcomp.eq.11) then
cM        cintbg   =cintbg   *Mso4/Msv
cM        cintbg05 =cintbg05 *Mso4/Msv
cM        cintbg125=cintbg125*Mso4/Msv
cM      elseif(kcomp.eq.13) then
      if(kcomp.eq.0) then
        do i=0,imax
          if(r(i).le.rbcn) then
            rhorbc(i)=rhobc
          else
            rhorbc(i)=rhobc*(rbcn/r(i))**(3.0-fracdim(i))
          endif
c          write(30,*) r(i), rhorbc(i), fracdim(i)
        enddo
        cintbg   =0.0
        cintbg05 =0.99*1.0e-3*(4.0*pi/3.0)*r(28)**3.0
     $                        *(rhorbc(28)*dndlrk(28))*d 
        cintbg125=0.03*1.0e-3*(4.0*pi/3.0)*r(31)**3.0
     $                        *(rhorbc(31)*dndlrk(31))*d 
        do i=0,imax
          dcintbg=1.0e-3*p2*r(i)**3.0*(rhorbc(i)*dndlrk(i))*d  
          cintbg=cintbg+dcintbg
          if(i.le.27) cintbg05 =cintbg05 +dcintbg
          if(i.ge.32) cintbg125=cintbg125+dcintbg
        enddo
      endif
c*****************************************
c      write(*,*) 'Cbc, Coc, Csu12, Csu3 og Cbg ='
c      write(*,*) cintbc,  cintoc, cintsc, cintsa, cintbg       
 
cM     dry volume fractions for sulfate, vsi, soot, vbci, oc, voci,
c     dry volume fractions for H2SO4+(NH4)2SO4, vsi, soot, vbci, oc, voci,
c     and background aerosol, vai. Note that vsi+vbci+voci+vai=1.  
      do i=1,imax 
       vtot=cbc(i)/rhobc+coc(i)/rhooc+csu12(i)/rhosv+csu3(i)/rhos
csoa     $      +cbg(i)/rhob 
     $      +cbg(i)/rhobg 
csoa       vcbg(i)=(cbg(i)/rhob)/vtot
       vcbg(i)=(cbg(i)/rhobg)/vtot
       vcbc(i)=(cbc(i)/rhobc)/vtot
       vcoc(i)=(coc(i)/rhooc)/vtot
       vcsu12(i)=(csu12(i)/rhosv)/vtot
       vcsu3(i)=(csu3(i)/rhos)/vtot
       vcsu(i)=vcsu12(i)+vcsu3(i) 
       vai(i)=vcbg(i)                     ! background (sulfate, OC, BC, SS or DU, or a mixture of two both for kcomp=1&4)
       vaisave(i)=vai(i)                  ! saved value for scaling w.r.t. internally mixed background
       vbci(i)=vcbc(i)                    ! non-background BC
       voci(i)=vcoc(i)                    ! non-background OC
       vsi(i)=vcsu(i)                     ! non-background sulfate
      enddo
csoa      if(kcomp.eq.4) then
c     Take into account the BC-fraction in the background OC&BC(Ait) mode, 
c     assuming voci=vbci=0 from coag:
csoa       factor=(1.0-fac)/(1.0-fac*(1.0-rhooc/rhobc))
csoa      do i=1,imax                                      
csoa       vai(i)=factor*vaisave(i)           ! background OM... (1-factor)*vaisave is background BC   ! Blir feil, tror jeg!!!!!
csoa       vbci(i)=vaisave(i)-vai(i)          ! background BC
csoa       voci(i)=1.e-40                     ! non-background OC
csoa       vsi(i)=1.0-vai(i)-vbci(i)-voci(i)  ! non-background sulfate
csoa       vsi(i)=max(0.0, vsi(i))
csoa      enddo
c       write(*,*) fac, factor
csoa      else
csoa       factor=1.0
c       vombg=0.0
csoa      endif  ! kcomp
c      do i=1,imax 
c        write(60,*) r(i), vsi(i)
c        write(61,*) r(i), vbci(i)
c        write(62,*) r(i), voci(i)
c        write(63,*) r(i), vai(i)
c      enddo

      return
      end  
