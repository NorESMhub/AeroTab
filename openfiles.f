      subroutine openfiles(kcomp,iopt,ib)

c **********************************************************************************
c     Opening files for use as input to NorESM (only for kcomp=0-10). nkcomp*.out 
c     (full modified size distributions) are presently not used in NorESM.
c  
c     Created by Alf Kirkevåg.
c **********************************************************************************

      integer kcomp, iopt, ib

      if(iopt.eq.1) then

       if(ib.eq.31) then  ! SW aerocom optics only
        if(kcomp.eq.1) then
          call system('mv aerocomk1.out aerocomk1_old.out')
          open(9500, file='aerocomk1.out')
        elseif(kcomp.eq.2) then
          call system('mv aerocomk2.out aerocomk2_old.out')
          open(9500, file='aerocomk2.out')
        elseif(kcomp.eq.3) then
          call system('mv aerocomk3.out aerocomk3_old.out')
          open(9500, file='aerocomk3.out')
        elseif(kcomp.eq.4) then
          call system('mv aerocomk4.out aerocomk4_old.out')
          open(9500, file='aerocomk4.out')
        elseif(kcomp.eq.5) then
          call system('mv aerocomk5.out aerocomk5_old.out')
          open(9500, file='aerocomk5.out')
        elseif(kcomp.eq.6) then
          call system('mv aerocomk6.out aerocomk6_old.out')
          open(9500, file='aerocomk6.out')
        elseif(kcomp.eq.7) then
          call system('mv aerocomk7.out aerocomk7_old.out')
          open(9500, file='aerocomk7.out')
        elseif(kcomp.eq.8) then
          call system('mv aerocomk8.out aerocomk8_old.out')
          open(9500, file='aerocomk8.out')
        elseif(kcomp.eq.9) then
          call system('mv aerocomk9.out aerocomk9_old.out')
          open(9500, file='aerocomk9.out')
        elseif(kcomp.eq.10) then
          call system('mv aerocomk10.out aerocomk10_old.out')
          open(9500, file='aerocomk10.out')
        elseif(kcomp.eq.0) then
          call system('mv aerocomk0.out aerocomk0_old.out')
          open(9500, file='aerocomk0.out')
        endif        
          if(kcomp.eq.1) then
          call system('mv aerodryk1.out aerodryk1_old.out')
          open(9600, file='aerodryk1.out')
        elseif(kcomp.eq.2) then
          call system('mv aerodryk2.out aerodryk2_old.out')
          open(9600, file='aerodryk2.out')
        elseif(kcomp.eq.3) then
          call system('mv aerodryk3.out aerodryk3_old.out')
          open(9600, file='aerodryk3.out')
        elseif(kcomp.eq.4) then
          call system('mv aerodryk4.out aerodryk4_old.out')
          open(9600, file='aerodryk4.out')
        elseif(kcomp.eq.5) then
          call system('mv aerodryk5.out aerodryk5_old.out')
          open(9600, file='aerodryk5.out')
        elseif(kcomp.eq.6) then
          call system('mv aerodryk6.out aerodryk6_old.out')
          open(9600, file='aerodryk6.out')
        elseif(kcomp.eq.7) then
          call system('mv aerodryk7.out aerodryk7_old.out')
          open(9600, file='aerodryk7.out')
        elseif(kcomp.eq.8) then
          call system('mv aerodryk8.out aerodryk8_old.out')
          open(9600, file='aerodryk8.out')
        elseif(kcomp.eq.9) then
          call system('mv aerodryk9.out aerodryk9_old.out')
          open(9600, file='aerodryk9.out')
        elseif(kcomp.eq.10) then
          call system('mv aerodryk10.out aerodryk10_old.out')
          open(9600, file='aerodryk10.out')
        elseif(kcomp.eq.0) then
          call system('mv aerodryk0.out aerodryk0_old.out')
          open(9600, file='aerodryk0.out')
        endif        
       endif ! ib=31

       if(ib.ne.19) then  ! SW CAM optics only
        if(kcomp.eq.1) then
          call system('mv kcomp1.out kcomp1_old.out')
          open(9000, file='kcomp1.out')
        elseif(kcomp.eq.2) then
          call system('mv kcomp2.out kcomp2_old.out')
          open(9000, file='kcomp2.out')
        elseif(kcomp.eq.3) then
          call system('mv kcomp3.out kcomp3_old.out')
          open(9000, file='kcomp3.out')
        elseif(kcomp.eq.4) then
          call system('mv kcomp4.out kcomp4_old.out')
          open(9000, file='kcomp4.out')
        elseif(kcomp.eq.5) then
          call system('mv kcomp5.out kcomp5_old.out')
          open(9000, file='kcomp5.out')
        elseif(kcomp.eq.6) then
          call system('mv kcomp6.out kcomp6_old.out')
          open(9000, file='kcomp6.out')
        elseif(kcomp.eq.7) then
          call system('mv kcomp7.out kcomp7_old.out')
          open(9000, file='kcomp7.out')
        elseif(kcomp.eq.8) then
          call system('mv kcomp8.out kcomp8_old.out')
          open(9000, file='kcomp8.out')
        elseif(kcomp.eq.9) then
          call system('mv kcomp9.out kcomp9_old.out')
          open(9000, file='kcomp9.out')
        elseif(kcomp.eq.10) then
          call system('mv kcomp10.out kcomp10_old.out')
          open(9000, file='kcomp10.out')
        elseif(kcomp.eq.0) then
          call system('mv kcomp0.out kcomp0_old.out')
          open(9000, file='kcomp0.out')
        endif        
       endif

        if(kcomp.eq.1) then
          call system('mv nkcomp1.out nkcomp1_old.out')
          open(9001, file='nkcomp1.out')
        elseif(kcomp.eq.2) then
          call system('mv nkcomp2.out nkcomp2_old.out')
          open(9001, file='nkcomp2.out')
        elseif(kcomp.eq.3) then
          call system('mv nkcomp3.out nkcomp3_old.out')
          open(9001, file='nkcomp3.out')
        elseif(kcomp.eq.4) then
          call system('mv nkcomp4.out nkcomp4_old.out')
          open(9001, file='nkcomp4.out')
        elseif(kcomp.eq.5) then
          call system('mv nkcomp5.out nkcomp5_old.out')
          open(9001, file='nkcomp5.out')
        elseif(kcomp.eq.6) then
          call system('mv nkcomp6.out nkcomp6_old.out')
          open(9001, file='nkcomp6.out')
        elseif(kcomp.eq.7) then
          call system('mv nkcomp7.out nkcomp7_old.out')
          open(9001, file='nkcomp7.out')
        elseif(kcomp.eq.8) then
          call system('mv nkcomp8.out nkcomp8_old.out')
          open(9001, file='nkcomp8.out')
        elseif(kcomp.eq.9) then
          call system('mv nkcomp9.out nkcomp9_old.out')
          open(9001, file='nkcomp9.out')
        elseif(kcomp.eq.10) then
          call system('mv nkcomp10.out nkcomp10_old.out')
          open(9001, file='nkcomp10.out')
        elseif(kcomp.eq.0) then
          call system('mv nkcomp0.out nkcomp0_old.out')
          open(9001, file='nkcomp0.out')
        endif        

       if(ib.eq.19) then  ! LW optics only
        if(kcomp.eq.1) then
          call system('mv lwkcomp1.out lwkcomp1_old.out')
          open(9009, file='lwkcomp1.out')
        elseif(kcomp.eq.2) then
          call system('mv lwkcomp2.out lwkcomp2_old.out')
          open(9009, file='lwkcomp2.out')
        elseif(kcomp.eq.3) then
          call system('mv lwkcomp3.out lwkcomp3_old.out')
          open(9009, file='lwkcomp3.out')
        elseif(kcomp.eq.4) then
          call system('mv lwkcomp4.out lwkcomp4_old.out')
          open(9009, file='lwkcomp4.out')
        elseif(kcomp.eq.5) then
          call system('mv lwkcomp5.out lwkcomp5_old.out')
          open(9009, file='lwkcomp5.out')
        elseif(kcomp.eq.6) then
          call system('mv lwkcomp6.out lwkcomp6_old.out')
          open(9009, file='lwkcomp6.out')
        elseif(kcomp.eq.7) then
          call system('mv lwkcomp7.out lwkcomp7_old.out')
          open(9009, file='lwkcomp7.out')
        elseif(kcomp.eq.8) then
          call system('mv lwkcomp8.out lwkcomp8_old.out')
          open(9009, file='lwkcomp8.out')
        elseif(kcomp.eq.9) then
          call system('mv lwkcomp9.out lwkcomp9_old.out')
          open(9009, file='lwkcomp9.out')
        elseif(kcomp.eq.10) then
          call system('mv lwkcomp10.out lwkcomp10_old.out')
          open(9009, file='lwkcomp10.out')
        elseif(kcomp.eq.0) then
          call system('mv lwkcomp0.out lwkcomp0_old.out')
          open(9009, file='lwkcomp0.out')
        endif        
       endif  ! ib=19

      else  ! iopt=0

        if(kcomp.eq.1) then
          call system('mv nkcomp1.out nkcomp1_old.out')
          open(9001, file='nkcomp1.out')
        elseif(kcomp.eq.2) then
          call system('mv nkcomp2.out nkcomp2_old.out')
          open(9001, file='nkcomp2.out')
        elseif(kcomp.eq.3) then
          call system('mv nkcomp3.out nkcomp3_old.out')
          open(9001, file='nkcomp3.out')
        elseif(kcomp.eq.4) then
          call system('mv nkcomp4.out nkcomp4_old.out')
          open(9001, file='nkcomp4.out')
        elseif(kcomp.eq.5) then
          call system('mv nkcomp5.out nkcomp5_old.out')
          open(9001, file='nkcomp5.out')
        elseif(kcomp.eq.6) then
          call system('mv nkcomp6.out nkcomp6_old.out')
          open(9001, file='nkcomp6.out')
        elseif(kcomp.eq.7) then
          call system('mv nkcomp7.out nkcomp7_old.out')
          open(9001, file='nkcomp7.out')
        elseif(kcomp.eq.8) then
          call system('mv nkcomp8.out nkcomp8_old.out')
          open(9001, file='nkcomp8.out')
        elseif(kcomp.eq.9) then
          call system('mv nkcomp9.out nkcomp9_old.out')
          open(9001, file='nkcomp9.out')
        elseif(kcomp.eq.10) then
          call system('mv nkcomp10.out nkcomp10_old.out')
          open(9001, file='nkcomp10.out')
        endif        

        if(kcomp.eq.1) then
          call system('mv logntilp1.out logntilp1_old.out')
          open(9003, file='logntilp1.out')
        elseif(kcomp.eq.2) then
          call system('mv logntilp2.out logntilp2_old.out')
          open(9003, file='logntilp2.out')
        elseif(kcomp.eq.3) then
          call system('mv logntilp3.out logntilp3_old.out')
          open(9003, file='logntilp3.out')
        elseif(kcomp.eq.4) then
          call system('mv logntilp4.out logntilp4_old.out')
          open(9003, file='logntilp4.out')
        elseif(kcomp.eq.5) then
          call system('mv logntilp5.out logntilp5_old.out')
          open(9003, file='logntilp5.out')
        elseif(kcomp.eq.6) then
          call system('mv logntilp6.out logntilp6_old.out')
          open(9003, file='logntilp6.out')
        elseif(kcomp.eq.7) then
          call system('mv logntilp7.out logntilp7_old.out')
          open(9003, file='logntilp7.out')
        elseif(kcomp.eq.8) then
          call system('mv logntilp8.out logntilp8_old.out')
          open(9003, file='logntilp8.out')
        elseif(kcomp.eq.9) then
          call system('mv logntilp9.out logntilp9_old.out')
          open(9003, file='logntilp9.out')
        elseif(kcomp.eq.10) then
          call system('mv logntilp10.out logntilp10_old.out')
          open(9003, file='logntilp10.out')
        endif        

      endif  ! iopt

      return
      end
