!--------------------------------------------------------------------------------
!
!  Copyright (C) 2017  L. J. Allen, H. G. Brown, A. J. D’Alfonso, S.D. Findlay, B. D. Forbes
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!   
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!                       
!--------------------------------------------------------------------------------

 module m_crystallography

     contains

      SUBROUTINE CRYST(A0,DEG,SS)
      use m_precision
      !
      !  This subroutine calculates the triclinic information and stores it
      !  in SS(7).
      !
      !  A0 contains a, b and c in Angs, DEG
      !  contains alpha, beta and gamma in degrees, and SS is an array to
      !  carry the triclinic information.
      !
      implicit none
      !implicit double precision (a-h,o-z)
      integer(4) i
      real(fp_kind) ang,derad,D,COMP
      real(fp_kind) A0(3), DEG(3), SS(7), C(3), S(3)
      DATA COMP / 1.0e-04_fp_kind /
      derad=atan(1.0_fp_kind)*4.0_fp_kind/180.0_fp_kind
      
      if(a0(1).le.0.0_fp_kind) go to 99
      do  i = 1, 3
            ang = deg(i) * derad
            c(i) = cos( ang )
            s(i) = sin( ang )

            if(abs(c(i)).lt.comp) then
                  c(i) = 0.0_fp_kind
                  s(i) = 1.0_fp_kind
        endif
      enddo

      SS(1)  =  (A0(2) * A0(3) * S(1)) ** 2.0_fp_kind
      SS(2)  =  (A0(1) * A0(3) * S(2)) ** 2.0_fp_kind
      SS(3)  =  (A0(1) * A0(2) * S(3)) ** 2.0_fp_kind
      D      =  A0(1) * A0(2) * A0(3)
      SS(4)  =  D * A0(3) * (C(1) * C(2) - C(3))
      SS(5)  =  D * A0(1) * (C(2) * C(3) - C(1))
      SS(6)  =  D * A0(2) * (C(3) * C(1) - C(2))
      SS(7)  =  D * sqrt( 1._fp_kind - C(1) ** 2.0_fp_kind - C(2) ** 2.0_fp_kind &
      & - C(3) ** 2.0_fp_kind + 2.0_fp_kind * C(1) * C(2) * C(3) )

99    continue
      return
      end
!--------------------------------------------------------------------------------------
	subroutine zone(ig1,ig2,izone)
      !
      ! Derived from zone axis law CJR 29/3/90
      ! Thanks to Peter Miller. Find [uvw] zone
      ! from input ig1(3). ig2(3) is rotated clockwise from ig1.
      !

      use m_precision
      implicit none
      integer(4) i,ic,izz
      integer(4) izone(3),ig1(3),ig2(3),max
      real(fp_kind) zzz,diff
      izone(1) = ig1(2) * ig2(3) - ig1(3) * ig2(2)
      izone(2) = ig1(3) * ig2(1) - ig1(1) * ig2(3)
      izone(3) = ig1(1) * ig2(2) - ig1(2) * ig2(1)
      max = 0
      do i = 1, 3
            izz = abs(izone(i))
	      max = max0(max,izz)
	enddo
      ic = 1
   21 continue
      do i = 1, 3
            zzz = float(ic * izone(i))
            diff = zzz / max - ic * izone(i) / max
            if(abs(diff).gt.0.01_fp_kind) then
		      ic = ic + 1
		      go to 21
            endif
	enddo
      do i = 1,3
		  izone(i) = izone(i) * float(ic)/float(max)
      enddo
	return
      end

!--------------------------------------------------------------------------------------
      subroutine subuvw(hkl,ruvw,a0,deg,ss)
      !
      ! This subroutine finds the real space vector [uvw] which is
      ! parallel to an input reciprocal lattice vector (hkl), and
      ! its magnitude is 1 / | hkl |.
      !
      use m_precision
      implicit none

      integer(4) i,j
      real(fp_kind) ss(7),ruvw(3),fuvw(3),a0(3),deg(3),auvw
      real(fp_kind) amaxim,amaxr,factor,ahkl
      integer(4) hkl(3)
      fuvw(1) = hkl(1)*ss(1)+hkl(2)*ss(4)+hkl(3)*ss(6)
      fuvw(2) = hkl(1)*ss(4)+hkl(2)*ss(2)+hkl(3)*ss(5)
      fuvw(3) = hkl(1)*ss(6)+hkl(2)*ss(5)+hkl(3)*ss(3)
      amaxim=0.0_fp_kind
      do i=1,3
		if(abs(fuvw(i)).gt.amaxim) then
			j = i
			amaxim = abs(fuvw(i))
			ruvw(j) = 1.0_fp_kind
		endif
	  enddo
      do i = 1,3
       ruvw(i) = ruvw(j)*fuvw(i)/fuvw(j)
      enddo
      ahkl = trimi(hkl,ss)
      auvw = rsd(ruvw,a0,deg)
      factor = auvw * ahkl

      amaxr = 0.0d0
      do i = 1,3
       ruvw(i) = ruvw(i)/factor
       if(abs(ruvw(i)).gt.amaxr) amaxr = abs(ruvw(i))
	  enddo

      do  i = 1,3
		if(abs(ruvw(i)/amaxr).lt.0.00001_fp_kind) ruvw(i) = 0.0_fp_kind
	  enddo
!	----------------------------------------------------------------
!	Bug fix inserted 4/6/99. LJA. Otherwise for 0 -2 2, for example,
!	the real space vector was in the wrong direction. If the first
!	nonzero component is positive then there was no problem.
!	----------------------------------------------------------------


	  if (fuvw(j).lt.0.0_fp_kind) then
		do i=1,3
			ruvw(i) = -ruvw(i)
		enddo
	  endif

      return
      end

!--------------------------------------------------------------------------------------
      subroutine subhkl(izone,gg,a0,deg,ss)
      !
      ! This subroutine finds the reciprocal space vector (gg) which is
      ! parallel to an input real space (integer) lattice vector [izone], and
      ! its magnitude is 1 / | izone |.
      !
      use m_precision
      implicit none

      real(fp_kind) gg(3),fuvw(3),a0(3),deg(3),ruvw(3),ss(7),c(3),auvw
      real(fp_kind) amaxim,radeg,ahkl,amaxg,factor,sign
      integer(4) izone(3)
      integer(4) i,j,ihit

      radeg=180.0_fp_kind/(atan(1.0_fp_kind)*4.0_fp_kind)  
      do i = 1,3
		c(i) = cos( deg(i) / radeg)
		if(c(i).gt.0.9999_fp_kind) c(i) = 1.0_fp_kind
	  enddo
      fuvw(1) = izone(1) * a0(1) ** 2.0_fp_kind +		  &
     &          izone(2) * a0(1) * a0(2) * c(3) + &
     &          izone(3) * a0(3) * a0(1) * c(2)	  
      fuvw(2) = izone(2) * a0(2) ** 2.0_fp_kind +		  &
     &          izone(3) * a0(2) * a0(3) * c(1) + &
     &          izone(1) * a0(1) * a0(2) * c(3)
      fuvw(3) = izone(3) * a0(3) ** 2.0_fp_kind +		  &
     &          izone(1) * a0(3) * a0(1) * c(2) + &
     &          izone(2) * a0(2) * a0(3) * c(1)

      amaxim=0.0_fp_kind
      ihit = 0

      do i = 1, 3
       ruvw(i) = dble(izone(i))
       if(abs(fuvw(i)).gt.amaxim) then
		j = i
		ihit = 1
		amaxim = abs(fuvw(i))
		gg(j) = 1.0_fp_kind
       endif
      enddo

      if(ihit.eq.1) then
       do i = 1,3
	gg(i) = gg(j) * fuvw(i) / fuvw(j)
	   enddo

       ahkl = trimr(gg,ss)
       auvw = rsd(ruvw,a0,deg)
       factor = auvw * ahkl

       sign = 1.0_fp_kind
       if(gg(j)/fuvw(j).lt.0.0_fp_kind) sign = -1.0_fp_kind
	       amaxg = 0.0_fp_kind
		   do i = 1,3
			gg(i) = sign * gg(i) / factor
			if(abs(gg(i)).gt.amaxg) amaxg = abs(gg(i))
		   enddo

		   do i = 1,3
			  if(abs(gg(i)/amaxg).lt.0.00001_fp_kind) gg(i) = 0.0_fp_kind
		   enddo
	   endif

	   if(ihit.eq.0) then
       do i = 1,3
			gg(i) = 0.0_fp_kind
	   enddo
      endif

      return
      end
!--------------------------------------------------------------------------------------
      subroutine rshkl(zone,gg,a0,deg,ss)
      !
      ! This subroutine finds the reciprocal space vector (gg) which is
      ! parallel to an input real space (real) lattice vector [zone], and
      ! its magnitude is 1 / | zone |.
      !
      use m_precision
      implicit none

      real(fp_kind) gg(3),fuvw(3),a0(3),deg(3),ruvw(3),ss(7),c(3)
      real(fp_kind) zone(3),auvw,ahkl,amaxg,amaxim,factor,radeg,sign
      integer(4) i,j

      radeg=180.0_fp_kind/(atan(1.0_fp_kind)*4.0_fp_kind)  

      do i = 1,3
		c(i) = cos( deg(i) / radeg)
		if(c(i).gt.0.99990_fp_kind) c(i) = 1.00_fp_kind
	  enddo
      fuvw(1) = zone(1) * a0(1) ** 2.00_fp_kind +		&
     &          zone(2) * a0(1) * a0(2) * c(3) +&
     &          zone(3) * a0(3) * a0(1) * c(2)
      fuvw(2) = zone(2) * a0(2) ** 2.00_fp_kind +		&
     &          zone(3) * a0(2) * a0(3) * c(1) +&
     &          zone(1) * a0(1) * a0(2) * c(3)
      fuvw(3) = zone(3) * a0(3) ** 2.00_fp_kind +		&
     &          zone(1) * a0(3) * a0(1) * c(2) +&
     &          zone(2) * a0(2) * a0(3) * c(1)

      amaxim=0.00_fp_kind

      do i=1,3
		ruvw(i) = zone(i)
		if(abs(fuvw(i)).gt.amaxim) then
			j = i
			amaxim = abs(fuvw(i))
			gg(j) = 1.00_fp_kind
		endif
	  enddo

      do i = 1,3
       gg(i) = gg(j) * fuvw(i) / fuvw(j)
      enddo

      ahkl = trimr(gg,ss)
      auvw = rsd(ruvw,a0,deg)
      factor = auvw * ahkl

      sign = 1.00_fp_kind
      if(gg(j)/fuvw(j).lt.0.00_fp_kind) sign = -1.00_fp_kind

      amaxg = 0.00_fp_kind
      do i = 1,3
       gg(i) = sign * gg(i) / factor
       if(abs(gg(i)).gt.amaxg) amaxg = abs(gg(i))
	  enddo

      do i = 1,3
		if(abs(gg(i)/amaxg).lt.0.000010_fp_kind) gg(i) = 0.00_fp_kind
	  enddo

      return
      end

!--------------------------------------------------------------------------------------
      subroutine angle(ig1,ig2,ss,thetad)
      !
      ! This subroutine finds the angle (thetad) in degrees between two input
      ! reciprocal lattice vectors g1 and g2. The angle lies between
      ! 0 and 180 degrees, and G2 is rotated by this angle clockwise
      ! from the G1 direction.
      use m_precision
      implicit none

      real(fp_kind) h12(3), h1(3), h2(3), ss(7)
      real(fp_kind) thetad,ag1,ag2,deg1,deg2
      real(fp_kind) pi
      integer(4) ig1(3), ig2(3)
      integer(4) i,j

      pi=atan(1.00_fp_kind)*4.00_fp_kind
      ag1 = trimi(ig1,ss)
      ag2 = trimi(ig2,ss)
      if(ag1.eq.0.0_fp_kind.or.ag2.eq.0.0_fp_kind) then
            write(6,101)
  101       format(' Error in Angle - one vector has zero magnitude')
            go to 99
      endif
      do i = 1,3
            h1(i) = float(ig1(i))/ag1
            h2(i) = float(ig2(i))/ag2
	      h12(i) = h1(i) + h2(i)
	enddo
      if(h12(1).eq.0.0_fp_kind.and.h12(2).eq.0.0_fp_kind.and.h12(3).eq.0.0_fp_kind) then
            thetad = 180.0_fp_kind
            go to 99
      endif
      deg1 = acos(cosanr(h1,h12,ss)) * 180.0_fp_kind / pi
      deg2 = acos(cosanr(h12,h2,ss)) * 180.0_fp_kind / pi
      thetad = deg1 + deg2
      if(abs(thetad-90.0_fp_kind).lt.0.0001_fp_kind) thetad = 90.0_fp_kind
   99 return
      end

!--------------------------------------------------------------------------------------      
      function trimr(A,SS)

      !  This function returns the magnitude of a NON INTEGER
      !  reciprocal lattice vector
      !  A(3), in A-1. SS(7) contains triclinic information.

      use m_precision
      implicit none

      real(fp_kind) A(3), SS(7), X1, X2, X3
      real(fp_kind) Y1, Y2, Y3, XX, YY
      real(fp_kind) trimr

      X1 = SS(1) * A(1) * A(1)
      X2 = SS(2) * A(2) * A(2)
      X3 = SS(3) * A(3) * A(3)
     
      Y1 = 2 * SS(4) * A(1) * A(2)
      Y2 = 2 * SS(5) * A(2) * A(3)
      Y3 = 2 * SS(6) * A(1) * A(3)
      XX = X1 + X2 + X3
      YY = Y1 + Y2 + Y3
      trimr = sqrt( abs( XX + YY ))/ SS(7)
      RETURN
      END

!--------------------------------------------------------------------------------------
      function trimi(A,SS)

      !  This function returns the magnitude of an INTEGER
      !  reciprocal lattice vector
      !  A(3), in A-1. SS(7) contains triclinic information.

      use m_precision
      implicit none

      real(fp_kind) SS(7), X1, X2, X3
      real(fp_kind) Y1, Y2, Y3, XX, YY
      real(fp_kind) trimi
      integer(4) A(3)

      X1 = SS(1) * dble(A(1)) * dble(A(1))
      X2 = SS(2) * dble(A(2)) * dble(A(2))
      X3 = SS(3) * dble(A(3)) * dble(A(3))

      Y1 = 2 * SS(4) * dble(A(1) * A(2))
      Y2 = 2 * SS(5) * dble(A(2) * A(3))
      Y3 = 2 * SS(6) * dble(A(1) * A(3))
      XX = X1 + X2 + X3
      YY = Y1 + Y2 + Y3
      trimi = sqrt( abs( XX + YY ))/ SS(7)
      RETURN
      END

!--------------------------------------------------------------------------------------
      FUNCTION RSD(Z,A0,DEG)

      !     This function returns the REAL SPACE DISTANCE of vector Z(3),
      !     Triclinic information is contained in A0 and DEG.

      use m_precision
      implicit none

      real(fp_kind) Z(3), A0(3), DEG(3), c(3)
      real(fp_kind) derad, one, two
      real(fp_kind) rsd
      integer(4) i

      derad=(atan(1.0_fp_kind)*4.0_fp_kind)/180.0_fp_kind
      do i = 1, 3
            c(i) = cos(derad * deg(i))
            if(c(i).gt.0.9999_fp_kind) c(i) = 1.0_fp_kind
      enddo
      

      ONE = 0.0_fp_kind
      do I = 1, 3
            ONE = (A0(I) * Z(I)) ** 2.0_fp_kind + ONE
	enddo

      TWO = 2.0_fp_kind * A0(2) * A0(3) * Z(2) * Z(3) * c(1)  + &
      &2.0_fp_kind * A0(3) * A0(1) * Z(3) * Z(1) * c(2)  + &
      &2.0_fp_kind * A0(1) * A0(2) * Z(1) * Z(2) * c(3)

      RSD = sqrt(ONE + TWO)
      RETURN
      END
!--------------------------------------------------------------------------------------
      FUNCTION COSANR(A,B,SS)

      !  This function returns the cosine of the angle between NON INTEGER
      !  reciprocal lattice vectors A(3) and B(3).
      !  SS(7) contains triclinic data.

      use m_precision
      implicit none

      real(fp_kind) A(3), B(3), SS(7)
      real(fp_kind) F, G, Y1, Y2, Y3, YY, Z1, Z2, Z3, ZZ
      real(fp_kind) prelim
      real(fp_kind) cosanr

      F = TRIMR( A, SS )
      G = TRIMR( B, SS )
      Y1 = SS(1) * A(1) * B(1)
      Y2 = SS(2) * A(2) * B(2)
      Y3 = SS(3) * A(3) * B(3)
      YY = Y1 + Y2 + Y3
      Z1 = SS(4) * (A(1) * B(2) + B(1) * A(2))
      Z2 = SS(5) * (A(2) * B(3) + B(2) * A(3))
      Z3 = SS(6) * (A(3) * B(1) + B(3) * A(1))
      ZZ = Z1 + Z2 + Z3
      prelim = (YY + ZZ) / (SS(7) ** 2.0_fp_kind * F * G)
      if(abs(prelim).gt.1.0_fp_kind) then
            cosanr = prelim/abs(prelim)
      else
            cosanr = prelim
      endif
      RETURN
      END


      !--------------------------------------------------------------------------------
      !
      ! This subroutine finds the angle (thetad) in degrees between two input
      ! real reciprocal lattice vectors g1 and g2. The angle lies between
      ! 0 and 180 degrees, and G2 is rotated by this angle clockwise
      ! from the G1 direction.

      subroutine angler(g1,g2,ss,thetad)
      use m_precision
      implicit none

      real(fp_kind) h12(3), h1(3), h2(3), ss(7)
      real(fp_kind) g1(3), g2(3)
      real(fp_kind) ag1, ag2, deg1, deg2, pi, thetad
      integer(4) i

      pi=4.0_fp_kind*atan(1.0_fp_kind)

	ag1 = trimr(g1,ss)
	ag2 = trimr(g2,ss)
	if(ag1.eq.0.0.or.ag2.eq.0.0) then
	write(6,101)
101   format(' Error in Angle - one vector has zero magnitude')
	go to 99
	endif
	do i = 1,3
	      h1(i) = g1(i)/ag1
	      h2(i) = g2(i)/ag2
            h12(i) = h1(i) + h2(i)
	enddo
      if(h12(1).eq.0.and.h12(2).eq.0.and.h12(3).eq.0) then
	thetad = 180.0_fp_kind
	go to 99
	endif
	deg1 = acos(cosanr(h1,h12,ss)) * 180.0_fp_kind / pi
	deg2 = acos(cosanr(h12,h2,ss)) * 180.0_fp_kind / pi
	thetad = deg1 + deg2
	if(abs(thetad-90).lt.0.0001_fp_kind) thetad = 90.0_fp_kind
 99   return
	end
	
end module m_crystallography

