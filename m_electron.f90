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

module m_electron
    
    implicit none
    
    contains
    
    
    
      function wavev(e)
      !  this function returns the wavevector in one over lambda, in a-1,
      !  for an input electron energy e in ev.
      
      use m_precision
      
      implicit none
      
      real(fp_kind) c1,c2,e
      real(fp_kind) wavev
      data c1, c2 / 9.78475598e-07_fp_kind, 12.2642596_fp_kind /
      
      wavev = sqrt( e + c1 *e ** 2.0_fp_kind ) / c2
      
      end function 
   
    
    
    
   end module
   