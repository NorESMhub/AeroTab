csoa      subroutine koehler (d, imax, r, vsi, vbci, voci, rh, f, fm, 
      subroutine koehler (d, imax, r, vsi, vbci, voci, vombg, vbcbg, 
     $           rh, f, fm, itot, faq, kcomp, iopt, supers, rsup,
     $           xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)

c **********************************************************************************
c     Created by Alf Kirkevåg.
c **********************************************************************************

c     solves the koehler equation to find wet/dry particle radii for
c     a given relative humidity, rh, or find the critical radius for 
c     CCN activation, rsup, given the supersaturation, supers
c     Note: minimum rsup obtainable for pure mineral aerosol is
c     approximately 0.08um for the current radius-resolution,
c     i.e. 10 times the resolution in the main program.

      implicit none

ccccc6ccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      INTEGER i, i1, imax, imaxx, j, jstep, irh, ia, itot, 
     $ jtest, irtest, iopt, kcomp, isup, koehlerplot, jxbound, jmaxf
      REAL rhow, sigm, Rg, rh, rhum, rhumg, T, ai, e, pi, 
     $ Mw, d, drdrh, drhdr0, drdr0, dfdr0, rad, vbcrad, frr0, x,
     $ rk(0:500), r(0:100),  vsi(0:100), vbci(0:100), voci(0:100),
     $ vsk(0:500), vbck(0:500), vock(0:500), f(-1:100), 
     $ fk(-1:500), fm(-1:100), fmk(-1:500), radm, scmax(500), 
     $ supers, rsup, rdry, faq
      REAL xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss
csoa
      REAL vombg, vbcbg
csoa

c      Matemathical constants
      DATA e, pi / 2.71828182845905e0, 3.141592654e0 /
c      Physical constants
      DATA  sigm, Rg / 7.6e-2, 8.3143 / 
c      Temperature in the koehler calculations
      DATA  T / 2.7315e2 /
c      Molecular weight and mass density of water 
      DATA  Mw, rhow / 1.8016e1, 1.0e3 /

c     define a 10 times finer r-grid     
      imaxx=10*imax
      do i=0,imaxx
        rk(i)=10.0**(0.1*d*(i-10)-3.0)
      enddo

      if(rh.lt.0.05.or.iopt.eq.0) then
c     we assume that the growth factor = 1 for RH < 5%,
c     and also for iopt=0 (it will then not be used) 
        do i=0,imax
          f(i)=1.0
          fm(i)=1.0
        enddo
      endif
      if(rh.lt.0.05) goto 99 

c     find interpolated dry vulume fractions for the new grid
      do i=1,imaxx
        i1=int(0.1*i)        
        if(i1.eq.0) then
          vsk(i)=vsi(1)    
          vbck(i)=vbci(1)  
          vock(i)=voci(1) 
        elseif(i1.ge.1.and.i1.lt.imax) then
          vsk(i)=vsi(i1)+(vsi(i1+1)-vsi(i1))*0.1*(i-10*i1)
          vbck(i)=vbci(i1)+(vbci(i1+1)-vbci(i1))*0.1*(i-10*i1)
          vock(i)=voci(i1)+(voci(i1+1)-voci(i1))*0.1*(i-10*i1)
        else
          vsk(i)=vsi(imax)    
          vbck(i)=vbci(imax)  
          vock(i)=voci(imax) 
        endif
ctest
c          vsk(i)=0.0
c          vbck(i)=0.0
c          vock(i)=0.0 
ctest
c          write(113,*) rk(i), vsk(i)
c          write(114,*) rk(i), vock(i)
c          write(115,*) rk(i), vbck(i)
      enddo

      if(iopt.eq.1) then ! for hygroscopic growth calculations      

c      calculate wet radii rad(i)=r(rh), such that f(rh)=rad(i)/rk(j),
c      i.e. given a dry radius rk, the wet radius rad is found by
c      solving the Koehler equation
       jstep=10
       do j=jstep,imaxx-jstep,jstep   ! dry radius index 
         i=j
         rhum=0.0
         jxbound=0
         do while(rhum.lt.rh.and.j.le.imaxx-jstep.and.i.le.imaxx-1)
           i=i+1                ! wet radius index 
           frr0=rk(i)/rk(j) 
c          mixsub calculates hygroscopic properties (given by x)
c          for an internally mixed aerosol
           call mixsub (frr0, itot, faq, Mw, rhow, 
csoa     $       j, vsk, vbck, vock, x, rh, kcomp,
     $       j, vsk, vbck, vock, vombg, vbcbg, x, rh, kcomp,
     $       xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)
           rhumg=rhum
c          the Koehler equation     
           rhum=e**(2e3*Mw*sigm/(Rg*T*rhow*rk(i))
     $       -x/((rk(i)/rk(j))**3-1.0))
           if(i.eq.imaxx.and.rhum.lt.rh) jxbound=jxbound+1
           if(jxbound.eq.1) jmaxf=j-jstep
c           if(j.eq.410) write(251,*) j, i, rhum    
         enddo
c        if rad < rk(imaxx) then use result from koehler formula 
         if(jxbound.eq.0) then  
           rad=((rk(i)-rk(i-1))*rh+rk(i-1)*rhum-rk(i)*rhumg)
     $         /(rhum-rhumg)       
           fk(j)=rad/rk(j)
c        but if rad > rk(imaxx), then use result for the largest
c        rad < rk(imaxx) in stead (a fair approximation) 
         else
           fk(j)=fk(jmaxf)
         endif
       enddo
       fk(imaxx)=fk(imaxx-jstep)

c      calculate also fm(i)=radm/rk(i), given a wet radius radm
       do i=jstep,imaxx,jstep   ! wet radius index 
         j=i
         rhum=0.0
         do while(rhum.lt.rh.and.j.ge.1)
           j=j-1                ! dry radius index 
           frr0=rk(i)/rk(j)
c          mixsub calculates hygroscopic properties (given by x)
c          for an internally mixed aerosol
           call mixsub (frr0, itot, faq, Mw, rhow, 
csoa     $       j, vsk, vbck, vock, x, rh, kcomp,
     $       j, vsk, vbck, vock, vombg, vbcbg, x, rh, kcomp,
     $       xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)
           rhumg=rhum
c          the Koehler equation     
           rhum=e**(2e3*Mw*sigm/(Rg*T*rhow*rk(i))
     $       -x/((rk(i)/rk(j))**3-1.0))    
         enddo
         radm=((rk(j+1)-rk(j))*rh+rk(j)*rhumg-rk(j+1)*rhum)
     $      /(rhumg-rhum)       
         fmk(i)=rk(i)/radm
       enddo

c      remap f and fm to the original resolution
       i=0
       do j=jstep,imaxx,jstep
         i=i+1
         f(i)=fk(j) 
         fm(i)=fmk(j) 
         if(i.gt.imax-2) f(i)=f(imax-2)
         if(i.gt.imax-2) fm(i)=fm(imax-2)
c         write(110,*) r(i), fm(i)
c         write(111,*) r(i), f(i)
       enddo

      elseif(iopt.eq.0) then

c          for plotting of the Koehler graph (only correct over hysteresis 
c          domain, since we use supers in stead of rh in the mixsub-call)
           koehlerplot=0
           if(koehlerplot.eq.1) then
c          chosen dry radius, rk(i)
           i=210
           rdry=rk(i)
           write(*,*) 'rdry =', rdry
            do j=1,imaxx,2
              if(rk(j).ge.rdry) then
              frr0=rk(j)/rdry
              call mixsub (frr0, itot, faq, Mw, rhow, 
csoa     $          i, vsk, vbck, vock, x, supers, kcomp,
     $          i, vsk, vbck, vock, vombg, vbcbg, x, supers, kcomp,
     $          xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)
c             the Koehler equation     
              rhum=e**(2e3*Mw*sigm/(Rg*T*rhow*rk(j))
     $          -x/((rk(j)/rdry)**3-1.0))    
              if(rhum.gt.1.0) rhum=100*(rhum-1.0)+1.0
c                write(112,*) rk(j), rhum
c                write(116,*) rk(j), x
              endif  
            enddo
           endif ! koehler-plot
c       Searching critical radius rsup for calculation of activated CCNs,
c       first by finding the top of the koehler curve for each rk(i) (dry radius).
c       Outer loop is over (increasing) dry radius, r(i)  
        do i=1,imaxx-1
          j=imaxx
c         initialization          
          rhumg=0.0
          rhum=0.0
          scmax(i)=0.0
          jtest=0
c          Inner loop is over (decreasing) wet radius, r(j)  
         do while(j.ge.i+2)
            j=j-1
c           for given particle size (index i) and growth factor frr0... 
            frr0=rk(j)/rk(i)
c           ...mixsub calculates hygroscopic properties (given by x) 
c           for an internally mixed aerosol...
            call mixsub (frr0, itot, faq, Mw, rhow, 
csoa     $        i, vsk, vbck, vock, x, supers, kcomp,
     $        i, vsk, vbck, vock, vombg, vbcbg, x, supers, kcomp,
     $        xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)
            rhumg=rhum
c           ... and the Koehler equation yield the corresponding 
c           equilibrium RH, rhum (= scmax at the top)      
            rhum=e**(2e3*Mw*sigm/(Rg*T*rhow*rk(j))
     $        -x/((rk(j)/rk(i))**3-1.0))
            if(rhum.gt.rhumg) then
              jtest=jtest+1
            else
              if(jtest.ge.1) then
                scmax(i)=rhumg
              else
                scmax(i)=1.0
              endif
              j=1
            endif   
          enddo
c          write(111,*) rk(i), scmax(i)
c          write(112,*) rk(i), rhum
c          write(116,*) rk(i), x
        enddo
c       scmax(i) gives RH for the top of the koehler curve given rk(i).
c       This last loop estimates the dry radius (rk) yielding scmax=supers
        irtest=0
        i=0
        do while(irtest.eq.0)
          i=i+1
          if(i.eq.imaxx+1) then
            write(*,*) 'no rsup found'
            stop
          endif
c          write(120,*) rk(i), scmax(i)
          if(scmax(i).le.supers.and.scmax(i).gt.1.0) then
            irtest=irtest+1
            if(irtest.eq.1) then
              if(i.gt.1) then
                rsup =
     $            (rk(i)*(supers-scmax(i-1))+rk(i-1)*(scmax(i)-supers))
     $           /(scmax(i)-scmax(i-1))       
              else
                rsup=rk(i)
              endif
            endif
          endif
        enddo

      endif   ! iopt=0

 99   return
      end  

