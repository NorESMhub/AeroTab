      subroutine mixsub (frr0, itot, faq, Mw, rhow,  
csoa     $ i, vsk, vbck, vock, x, rh, kcomp, 
     $ i, vsk, vbck, vock, vombg, vbcbg, x, rh, kcomp, 
     $ xbc, xdst, xoc, xs, xa, xss, rhda, rhca, rhdss, rhcss)
     
c **********************************************************************************
c     Created by Alf Kirkevåg.
c **********************************************************************************

c     Mixsub calculates hygroscopic properties (given by x) for an internal mixture 
c     of aerosol components with different hygroscopicity (x*).
csoa  Internal mixture of two constituents in the background aerosol for kcomp = 4
csoa  and new kcomp = 1 (including SOA nucleation) is now taken into account.

      implicit none

      INTEGER itot, i, kcomp
      REAL rhow, rhosl, ai, Mw, Ms, frh, frr0, x, xws, xss, xbc, xoc,
     $ xs, xa, xm, rh, vsk(0:10000), vbck(0:10000), vock(0:10000), 
     $ xwaso, xdst, xbg, faq
      REAL rhda, rhca, rhdss, rhcss
csoa
      REAL vombg, vbcbg
csoa

c     Set the hygroscopicity of the background aerosol 
csoa      if(kcomp.eq.1.or.kcomp.eq.5) then
      if(kcomp.eq.1) then
c       Sulfuric acid
csoa    internally mixed with OM (as SOA) (using volume mixing approximation):
        xbg=xs*(1-vombg)+xoc*vombg
ct         xbg=xs   ! this is the correct value for H2SO4
ct         xbg=xa   ! this gives hygroscopicity for ammonium sulfate instead of H2SO4 (just for testing & plotting purposes)
      elseif(kcomp.eq.2) then 
c        Hydrophobic BC:
         xbg=xbc
      elseif(kcomp.eq.3) then
c        Organic Carbon:                      
         xbg=xoc
      elseif(kcomp.eq.4) then
c        Organic Carbon:                      
csoa         xbg=xoc
         xbg=xoc*(1-vbcbg)+xbc*vbcbg
csoa
      elseif(kcomp.eq.5) then 
c        Sulfuric acid:
         xbg=xs   ! this is the correct value for H2SO4
ct         xbg=xa   ! this gives hygroscopicity for ammonium sulfate instead of H2SO4 (just for testing & plotting purposes)
      elseif(kcomp.eq.6.or.kcomp.eq.7) then 
c        Mineral dust:
         xbg=xdst
      elseif(kcomp.eq.8.or.kcomp.eq.9.or.kcomp.eq.10) then 
c        Sea-salt:
         xbg=xss
      endif 

c     Hygroscopicity x for internally mixed aerosol is calculated by using the 
c     volume mixing assumprion, based on background aerosol, sulfate, soot and oc.
c     note: internally mixed sulphate is here assumed to be ammoniumsulfate,
c     except for mode 1-4, where all SO4 is H2SO4 instead.
      if(kcomp.ge.1.and.kcomp.le.10) then
       if(itot.eq.0) then
         x=xbg
       else    ! internal mixture                                                     Prøv med coating-antakelser her!!!!!
         x=(1.0-vsk(i)-vbck(i)-vock(i))*xbg
     $    +vsk(i)*(faq*xa+(1.0-faq)*xs)+vbck(i)*xbc+vock(i)*xoc 
       endif 
c     only sulfate or soot:
      else 
        write(*,*) 'kcomp = 1-10 only'
        stop
      endif

c      write(117,*), i, x, rh

      return
      end  

