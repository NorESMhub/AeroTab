      subroutine tabrefind (kcomp, ib, xlam, cref) 

c **********************************************************************************
c     Created by Alf Kirkevåg.
c **********************************************************************************

c      Here wavelength dependent complex rafractive indices (cref) 
c      are found from tabulated values for each aerosol component. 

      implicit none

      INTEGER  j, kcomp, ib, iband
      REAL     xlam(31), lam, lamg, nr, nrg, ni, nig
      COMPLEX  crin, cref(5,31)

c     initializing cref-array:
        do iband = 1, 31
      do j=1,5
        cref(j,iband) = 0.0
      enddo
        enddo 

        do iband = 1, ib

      do 300 j=1,5

c      refractive indices for background aerosol 
csoa   Here bacground aerosol does not take into account internal mixing,
csoa   so this must be taken into account separately in refind.f for modes
csoa   which have internal mixtures of two constituents! 
       if(j.eq.1) then
         if(kcomp.eq.1.or.kcomp.eq.5) then
           open(10, file='input/suso_gads.inp', status='old') 
         elseif(kcomp.eq.2) then
           open(10, file='input/sot_janzen.inp', status='old') 
         elseif(kcomp.eq.3.or.kcomp.eq.4) then
           open(10, file='input/waso_gads.inp', status='old') 
         elseif(kcomp.eq.6.or.kcomp.eq.7) then
!orig           open(10, file='input/mineral_mix.inp', status='old') 
           open(10, file='input/mineral_gads.inp', status='old') 
!test
         elseif(kcomp.eq.8.or.kcomp.eq.9.or.kcomp.eq.10) then
           open(10, file='input/ss_gads.inp', status='old')
         else
           goto 300
         endif
c      refractive indices for components to be internally mixed 
c      with the background aerosol
       elseif(j.eq.2) then
         open(10, file='input/suso_gads.inp', status='old')
       elseif(j.eq.3) then
         open(10, file='input/sot_janzen.inp', status='old')
       elseif(j.eq.4) then
         open(10, file='input/water_gads.inp', status='old')
       elseif(j.eq.5) then ! assumed OC refractive index
         open(10, file='input/waso_gads.inp', status='old')
       endif

       lamg=0
       nrg=0
       nig=0

       read(10,*)
       read(10,*)
       read(10,*)
       read(10,*)
       read(10,*)
 100   read(10,*) lam, nr, ni
       if(lam.gt.xlam(iband)) then
         cref(j,iband)=
     $  (1e0,0e0)*(xlam(iband)*(nr-nrg)+(nrg*lam-nr*lamg))/(lam-lamg)
     $ -(0e0,1e0)*(xlam(iband)*(ni-nig)+(nig*lam-ni*lamg))/(lam-lamg)
         goto 200
       else           
         lamg=lam
         nrg=nr
         nig=ni
         goto 100
       endif
 200   close(10)
 
 300  continue

        enddo ! iband

	end










