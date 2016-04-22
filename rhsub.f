      subroutine rhsub (imax, rh, d, r, rp, dndlrkny, vsi, vbci, voci,
csoa     $ vssol, vbcsol, vocsol, vasol, vw, fki, itot, rhos, rhosv, 
     $ fombg, fbcbg, vombg, vbcbg, 
     $ vssol, vbcsol, vocsol, vasol, vw, fki, itot, rhos, rhosv, 
     $ rhobc, rhooc, rhob, rhow, Ctot, kcomp, iopt, supers, rsup, 
csoa     $ ismolarh, iSOA, cat, fac, fabc, faq, CCN,
     $ ismolarh, cat, fac, fabc, faq, CCN,
     $ xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)

c **********************************************************************************
c     Created by Alf Kirkevåg.
c **********************************************************************************

c     Hygroscopic growth is taken into account, either for the given 
c     relative humidity (if iopt=1), or for the given supersaturation 
c     (if iopt=0). New number and mass concentrations and volume fractions
c     are calculated, as well as numbers of activated CCN (if iopt=0). 

      implicit none

      INTEGER i, imax, j, jmax, itot, kcomp, iopt, ismolarh
      REAL dninc(0:100), dip(0:100), dndlrkny(0:100), dndlrccn(0:100),
     $ r(0:100), rp(0:100), rh, d, vbci(0:100), voci(0:100), 
     $ vsi(0:100), vssol(0:100), vbcsol(0:100), vocsol(0:100),
     $ vasol(0:100),vw(0:100), 
     $ fki(-1:100), fmax, rny, f(-1:100), fm(-1:100), vssolub(100), 
     $ vbcsolub(100), vocsolub(100), vasolub(100), dncny(0:100), Ctot, 
     $ dCtot, rhos, rhosv, rhobc, rhooc, rhob, rhow, rsup, supers,
     $ Nbak, CCN, wccn, cat, fac, fabc, faq, pi, e
      REAL xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss
csoa
      REAL fombg, fbcbg, vombg, vbcbg
csoa 
      PARAMETER  (pi=3.141592654, e=2.718281828)
csoa      LOGICAL khyd
cSOA
csoa      INTEGER iSOA
cSOA

csoa  "hydrophobic" (less than the added mass onto the) background components 
c     ( i.e. components that may have smaller f values for large r, despite 
c     potentially large f at smaller sizes due to internally mixed sulfate onto 
c     the smallest particles)
csoa      khyd=kcomp.ne.1.and.kcomp.ne.5 
     
c     save dndlrkny of the dry aerosol as dndlrccn, 
c     for use in the CCN calculations 
      do i=1,imax       
        dndlrccn(i)=dndlrkny(i)
      enddo                                

c     initialize wet volume fractions for sulfate, vssol, soot, vbcsol,
c     and background aerosol, vasol. 
csoa  Note that the background aerosol consists of an internal mixture
csoa  of two constituents for kcomp=1 (Sulfate and OM) and kcomp=4 (OM and BC)
      do i=1,imax
        vssol(i)=vsi(i)                      ! non-background sulfate
        vbcsol(i)=vbci(i)                    ! non-background BC
        vocsol(i)=voci(i)                    ! non-background OC
csoa        vasol(i)=1.0-vsi(i)-vbci(i)-voci(i)  ! background (sulfate, OC, BC, SS or DU, or a mixture of two if kcomp=1 or 4)
        vasol(i)=max(1.0-vsi(i)-vbci(i)-voci(i),0.0)  ! background (sulfate, OC, BC, SS or DU, or a mixture of two if kcomp=1 or 4)
      enddo

c     subroutine koehler solves the koehler equation to find wet particle 
c     radii for a given relative humidity, rh, or the critical radius for 
c     CCN activation, rsup, given the supersaturation, supers
csoa      call koehler (d, imax, r, vsi, vbci, voci, rh, f, fm, 
      call koehler (d, imax, r, vsi, vbci, voci, vombg, vbcbg, 
     $ rh, f, fm, itot, faq, kcomp, iopt, supers, rsup,
     $ xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)

csoa      if(khyd) fmax=1.0                         ! initialized fmax for a mode with "hydrophobic" background
      do i=-1,imax
        if(i.le.0) then
          f(i)=1.0
          fm(i)=1.0
        endif
c       fki is for use in sub-routine refind
        fki(i)=f(i)
csoa        if(khyd.and.fm(i).ge.fmax) fmax=fm(i)  ! calculated fmax value for a mode with "hydrophobic" background 
c        if(i.eq.12) write(90,*) rh, fki(i) ! r=0.0118 (ca)  mode 1 & 2
c        if(i.eq.14) write(90,*) rh, fki(i) ! r=0.022 (ca)   mode 8
c        if(i.eq.17) write(90,*) rh, fki(i) ! r=0.04 (ca)    mode 4
c        if(i.eq.20) write(90,*) rh, fki(i) ! r=0.075 (ca)   mode 5
c        if(i.eq.21) write(90,*) rh, fki(i) ! r=0.1 (ca) 
c        if(i.eq.22) write(90,*) rh, fki(i) ! r=0.13 (ca)    mode 9
c        if(i.eq.24) write(90,*) rh, fki(i) ! r=0.22 (ca)    mode 6
c        if(i.eq.29) write(90,*) rh, fki(i) ! r=0.63 (ca)    mode 7
c        if(i.eq.30) write(90,*) rh, fki(i) ! r=0.74 (ca)    mode 10
c        if(i.eq.43) write(90,*) rh, fki(i)
c        if(i.eq.30) write(91,*) rh, fm(i)  ! test 4nov2014 (ikke bruk denne)
      enddo
csoa      if(.not.khyd) fmax=f(imax)               ! calculated fmax value for a mode with "hygrophilic" background 
c      write(88,*) rh, fmax
c      write(90,*) rh, f(10)
c      write(94,*) rh, f(44)
c      write(95,*) rh, fm(10)
c      write(99,*) rh, fm(44)
c      write(*,*) rh, fmax
csoa
      fmax=1.0
      do i=1,imax
        if(fm(i).gt.fmax) then
          fmax=fm(i)
        endif
      enddo
csoa

c     the total iteration number jmax must be sufficiently large to 
c     satisfy the stability criterium for the continuity equation.
      jmax=int(log10(fmax)/d)+1
ctest
c     Note: when jmax is chosen large (larger than necessary), the solution 
c     becomes very diffusive
ctest
c      write(*,*) 'fmax, jmax =', fmax, jmax
      
c     determine the increment of log(r/um) at r=rp, i.e. in the center of 
c     the size bin, dip (chosen to be the same for every time step). 
      dip(0)=0.0
      do i=1,imax
        if(i.eq.imax) then
          dip(i)=log10(fm(i))/jmax
        else
          dip(i)=log10(0.5*(fm(i)+fm(i+1)))/jmax
        endif
      enddo

c     solve the continuity equations (with a simple upwind advection scheme,
c     or with corrective anti-diffusive steps from the Smolarkiewicz scheme) 
c     using jmax time steps/iterations for the size distribution, dndlrkny. 
c     Process specific wet volume fractions and mass concentrations are 
c     determined directly from the growth factor and the dry values

         do j=1,jmax

       do i=1,imax
         dninc(i)=-(dndlrkny(i)*dip(i)-dndlrkny(i-1)*dip(i-1))/d
       enddo
       do i=1,imax       
         dndlrkny(i)=dndlrkny(i)+dninc(i)
       enddo

       do i=1,imax       
          rny=r(i)*fm(i)**(-1.0/real(jmax))
         if(i.eq.1) then
           vssolub(i)=vssol(i)
           vbcsolub(i)=vbcsol(i)
           vocsolub(i)=vocsol(i)
           vasolub(i)=vasol(i)
         else
           vssolub(i)=(vssol(i-1)*log10(r(i)/rny)       
     $     +vssol(i)*log10(rny/r(i-1)))/d
           vbcsolub(i)=(vbcsol(i-1)*log10(r(i)/rny)       
     $     +vbcsol(i)*log10(rny/r(i-1)))/d
           vocsolub(i)=(vocsol(i-1)*log10(r(i)/rny)       
     $     +vocsol(i)*log10(rny/r(i-1)))/d
           vasolub(i)=(vasol(i-1)*log10(r(i)/rny)       
     $     +vasol(i)*log10(rny/r(i-1)))/d
         endif
       enddo
       do i=1,imax       
         vssol(i)=vssolub(i)*fm(i)**(-3.0/real(jmax))
         vbcsol(i)=vbcsolub(i)*fm(i)**(-3.0/real(jmax))
         vocsol(i)=vocsolub(i)*fm(i)**(-3.0/real(jmax))
         vasol(i)=vasolub(i)*fm(i)**(-3.0/real(jmax))
         vw(i)=1.0-min(vssol(i)+vbcsol(i)+vocsol(i)+vasol(i),1.0)
       enddo

      if(ismolarh.gt.0) then
c       Smolarkiewicz-scheme with ismolar corrective steps
        do i=1,imax
          dncny(i)=dndlrkny(i)
        enddo
        call smolar (ismolarh, imax, d, dncny, dip)  
        do i=1,imax
corig          dndlrkny(i)=dncny(i)
cfix  needed to avoid NaN with the new compiler on Precise
          dndlrkny(i)=max(1.e-80,dncny(i)) !test
        enddo
      endif

         enddo                     ! j-loop

c     volume fractions for sulfate, vssol, soot, vbcsol, oc, vocsol, 
c     background aerosol, vasol, and water, vw, after hygroscopic growth. 
c     Here vssol+vbcsol+vocsol+vasol+vw=1.
c      do i=1,imax
c        write(132,100) r(i), vssol(i)
c        write(133,100) r(i), vbcsol(i)
c        write(134,100) r(i), vocsol(i)
c        write(135,100) r(i), vasol(i)
c        write(136,100) r(i), vw(i)
c      enddo

c     condensed water contribution, dCtot, to the total aerosol 
c     concentration, Ctot (ug/m**-3) 
      do i=1,imax 
        dCtot=1.0e-3*(4.0*pi/3.0)*r(i)**3.0
     $        *(rhow*vw(i)*dndlrkny(i))*d
        Ctot=Ctot+dCtot                  
      enddo
c      write(*,*) 'wet Ctot =', Ctot
 
      CCN=0.0
      if(iopt.eq.0) then  
cNB     This part of the code is not updated: do not use without cleaning up and checking first! 
cNB:    It is Wrong for kcomp=8-10 at least...
c       calculation of activated CCN, i.e. the number of (dry) aerosol particles with r > rsup 
c        write(*,*) 'S, rsup =', 100*supers-100.0, rsup      
c        write(101,*) 100*supers-100.0, rsup      
        Nbak=0.0
        do i=1,imax       
          if(rp(i).ge.rsup) then
            if(rp(i-1).ge.rsup) then
              wccn=1.0
            else
              wccn=log10(rp(i)/rsup)/d
            endif
          else
            wccn=0.0
          endif
          CCN=CCN+wccn*dndlrccn(i)*d
          Nbak=Nbak+dndlrccn(i)*d          
c          write(35,*) r(i), dndlrccn(i), CCN
        enddo
c        write(*,*)
c        write(*,600) 'S, rsup, CCN, Nbak =', 
c     $   100*supers-100.0, rsup, CCN, Nbak
        write(*,700) 'S,rs,cat,fac,fabc,faq,CCN =', 
     $   100*supers-100.0,rsup,cat,fac,fabc,faq,CCN
c        write(*,*)
c       write the (tabulated) result to file 
ccccc6ccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
csoa        if(kcomp.ge.1.and.kcomp.le.3) then
        if(kcomp.eq.1) then
cSOA
csoa          if(kcomp.eq.1.and.iSOA.eq.1) then
           write(9002,500) kcomp, supers, fombg, cat, fac, CCN
        elseif(kcomp.eq.2.or.kcomp.eq.3) then
cSOA
csoa          else ! kcomp.eq.2.or.kcomp.eq.3
csoa           write(9002,400) kcomp, supers, cat, CCN
csoa          endif
           write(9002,410) kcomp, supers, cat, fac, CCN
        elseif(kcomp.eq.4) then
csoa          write(9002,500) kcomp, supers, cat, fac, faq, CCN
          write(9002,200) kcomp, supers, fbcbg, cat, fac, faq, CCN
        elseif(kcomp.ge.5.and.kcomp.le.10) then
          write(9002,200) kcomp, supers, cat, fac, fabc, faq, CCN
        else
          write(*,*) 'No valid kcomp>10, and 0 is assumed hydrofobic'
          stop
        endif

      endif
ccccc6ccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

 100  format(2(x,e10.4))
 200  format(I3,x,e13.6,6(x,e12.5))
 400  format(I3,x,e13.6,2(x,e12.5))
 410  format(I3,x,e13.6,3(x,e12.5))
 500  format(I3,x,e13.6,4(x,e12.5))
 600  format(A20,f9.2,x,f10.3,2(2x,e13.6))
 700  format(A27,f7.2,f8.3,x,e11.4,3f7.2,x,e10.3)

      return
      end

