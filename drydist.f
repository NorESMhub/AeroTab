      subroutine drydist(kcomp, Nnat, imini, imax, d, r, rk,  
     $ logsk, logs0, rhob, bcint, pi, dndlrk, ntot, Ctot)

c **********************************************************************************
c     Created by Alf Kirkevåg.
c **********************************************************************************

c     Dry background mode size distribution, dndlrk, and dry background
c     aerosol contribution to the mass concentration Ctot (ug/m**3). 
csoa  Note: this calculation of Ctot does not take into account that the 
csoa  backrgound aerosol can be an internal mixture of two constituents.
csoa  This must be compensated for (by scaling) other places in the code. 

      implicit none

      INTEGER kcomp, i, imini, imax
      REAL d, rk, logsk, nk, ntot, Nnat, rhob, bcint, pi, Ctot, dCtot 
      REAL logs0
      REAL r(0:100), dndlrk(0:100)

      ntot=0.0
      do i=0,imax
        nk=(1.0/logsk)*exp(-0.5*(log10(r(i)/rk)/logsk)**2.0)
        dndlrk(i)=Nnat*nk/sqrt(2.0*pi)
        ntot=ntot+dndlrk(i)*d
ctest        write(*,*) 'r, ntot =', i, r(i), ntot
c        write(14,*) r(i), dndlrk(i)
      enddo

      Ctot=0.0
      if(kcomp.ne.0) then
        do i=imini,imax  
ctst         if(r(i).ge.1.25) then  ! amount larger than 2.5um in diameter 
          dCtot=1.0e-3*(4.0*pi/3.0)*r(i)**3.0*(rhob*dndlrk(i))*d  
          Ctot=Ctot+dCtot
c          write(*,*) 'r, Ctot =', r(i), Ctot
ctst         endif
        enddo
      else
        dCtot=Nnat/(3e3*logs0/(sqrt(8.0*pi)*bcint))
        Ctot=Ctot+dCtot
      endif
c      write(*,*) 'Background mass conc. (normalized mode) =', Ctot
 
      return
      end 
