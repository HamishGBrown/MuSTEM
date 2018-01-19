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

module m_slicing
    use m_precision
    
    implicit none
    
    integer(4) :: n_slices                             !number of subslices in the unit cell
    integer(4) :: maxnat_slice                            !maximum number of atoms for a particular type in the supercell
    real(fp_kind), allocatable :: a0_slice(:,:)
    integer(4), allocatable :: nat_slice(:,:)             !number of atoms for each type in the slice (supercell)
    integer(4), allocatable :: nat_slice_unitcell(:,:)    !number of atoms for each type in the slice (unit cell)
    real(fp_kind), allocatable :: tau_slice(:,:,:,:)      !the atom positions for each potential subslice (supercell)
    real(fp_kind), allocatable :: tau_slice_unitcell(:,:,:,:) !the atom positions for each potential subslice (unit cell)
    real(fp_kind), allocatable :: prop_distance(:)        !Unit cell subslice propagation distance
    real(fp_kind), allocatable :: depths(:)               !Unit cell subslice depths (i.e. same as above)
    real(fp_kind), allocatable :: ss_slice(:,:)
    
    contains

    subroutine setup_slicing_depths()
    
        use m_user_input, only: get_input
        use m_precision, only: fp_kind

        implicit none
    
        integer :: i_slice
        
        write(*,*) '|-------------------------------|'
	    write(*,*) '|      Unit cell slicing        |'
	    write(*,*) '|-------------------------------|'
        write(*,*)
      
    22  write(6,23) char(143)
    23  format(' Do you wish to slice the unit cell in the beam direction?', /, &
              &' This may be a good idea if the unit cell is larger than 2 ', a1, /, &
              &' in the z-direction.', /, &
              &' <1> Yes <2> No ')
        call get_input('Slice unit cell <1> yes <2> no', i_slice)
        write(*,*) 
    
        if (i_slice.eq.1) then
            call calculate_depths_slicing
        
        elseif (i_slice.eq.2) then
            call calculate_depths_no_slicing
        
        else
            goto 22
        
        endif
    
    end subroutine
    
    subroutine calculate_depths_no_slicing
        implicit none
        
        n_slices = 1
        
        allocate(depths(2))
        
        depths(1) = 0.0_fp_kind
        depths(2) = 1.0_fp_kind
        
    end subroutine
    
    subroutine calculate_depths_slicing
        use m_user_input
        
        implicit none
        
        integer(4) i, ichoice

      5 write(*,*) 'Enter the number of slices per unit cell in the beam direction:'
	    call get_input("Number of distinct potentials ", n_slices)
        write(*,*)
        
        if (n_slices.eq.1) then
            call calculate_depths_no_slicing
                        
        elseif (n_slices.gt.1) then
            allocate(depths(n_slices+1))

            depths(n_slices+1) = 1.0_fp_kind
            
            write(6,10)
         10 format( ' You will now be asked to enter the depths at which slicing',/,&
                    &' will take place. These should be entered as fractions of',/,&
                    &' the unit cell. It is the front of the slice which should',/,&
                    &' be entered. (e.g. three even slices would be entered as',/,&
                    &' 0.0, 0.3333, 0.6666 ).', /)
	    
                    
    15      write(6,16) 
    16      format(' How would you like to specify the slice depths?', /, &
                  &' <1> Manually ', /, &
                  &' <2> Automatically (uniformly spaced)')
	        call get_input("<1> manual or <2> auto slicing", ichoice)
            write(*,*)
        
            if(ichoice.eq.1) then
                ! Manual slicing
            
	            do i = 1, n_slices
	                write(6,20,ADVANCE='NO') achar(13), i
        20          format( a1,' Enter the fractional depth of slice number ', i4, ':  ')
	                call get_input("depths", depths(i))
                    
	            enddo
                
            elseif (ichoice.eq.2) then
                ! Automatic slicing
            
    21          format( ' Fractional depth of slice ', i4, ':  ', f7.4)
                do i = 1, n_slices
                    depths(i) = float( (i-1)) / float(n_slices)
                    write(6,21) i, depths(i)
                enddo
                
            else
                ! Invalid slicing choice
                goto 15
                
            endif
            
        else
            ! Invalid number of slices
            goto 5
            
        endif
        
        write(*,*)
        
    end subroutine
    
    subroutine calculate_slices
        use global_variables, only: nat, ifactory, ifactorx, nt, a0, deg, tau, nm
        use m_crystallography, only: cryst
        
        implicit none
        
        integer :: j
        
        maxnat_slice = maxval(nat) * ifactory * ifactorx
        
        allocate(nat_slice(nt,n_slices))
        allocate(tau_slice(3,nt,maxnat_slice,n_slices))
        allocate(a0_slice(3,n_slices))
        allocate(ss_slice(7,n_slices))
        allocate(prop_distance(n_slices))
        allocate(nat_slice_unitcell(nt,n_slices))
        allocate(tau_slice_unitcell(3,nt,nm,n_slices))
        
        do j = 1, n_slices
          a0_slice(1,j) = a0(1)*ifactorx	
          a0_slice(2,j) = a0(2)*ifactory		
          
          if (j .eq. n_slices) then				
                a0_slice(3,j) = (1.0_fp_kind-depths(j))*a0(3)
          else
                a0_slice(3,j) = (depths(j+1)-depths(j))*a0(3)
          endif
          
          call cryst(a0_slice(:,j),deg,ss_slice(:,j))
          
          ! Make supercell cell subsliced
          call calculate_tau_nat_slice( depths(j),depths(j+1),tau,tau_slice(:,:,:,j),nm,nt,nat,nat_slice(:,j),ifactorx,ifactory )
	                
          ! Make unit cell subsliced
          call make_mod_tau_unitcell( depths(j),depths(j+1),tau,tau_slice_unitcell(:,:,:,j),nm,nt,nat,nat_slice_unitcell(:,j) )
        enddo
        
        prop_distance = a0_slice(3,:)
        
    end subroutine
    
	subroutine calculate_tau_nat_slice( depth1, depth2, tau, tau_slice, nm, nt, nat, nat_slice, ifactorx, ifactory )
   
    implicit none

	integer(4) nm, nt, ifactorx, ifactory
	integer(4) nat(nt)
	integer(4) nat_slice(nt)
	real(fp_kind) depth2, depth1, tau(3,nt,nm)
	real(fp_kind) tau_slice(3,nt,nm*ifactorx*ifactory)
	integer(4) i,j,jj, mm, nn
	real(fp_kind) diff, tol
	real(fp_kind) factorx, factory

	tol = 1.0d-4

	diff = depth2 - depth1

	factorx = float( ifactorx )
	factory = float( ifactory )

	do i = 1, nt
    
	   jj = 1
       
	   do mm = 1, ifactory
	   do nn = 1, ifactorx
	      do j = 1, nat(i)
          
	         if( (tau(3,i,j) .lt. (depth2-tol)) .and.(tau(3,i,j) .ge. (depth1-tol)) ) then
	            tau_slice(1,i,jj) = tau(1,i,j)/factorx + float(nn-1)/factorx
	            tau_slice(2,i,jj) = tau(2,i,j)/factory + float(mm-1)/factory
			    tau_slice(3,i,jj) = (tau(3,i,j)-depth1)/diff
	            jj = jj + 1
                
	         elseif(diff.eq.1.0) then
	            tau_slice(1,i,jj) = tau(1,i,j)/factorx + float(nn-1)/factorx
	            tau_slice(2,i,jj) = tau(2,i,j)/factory + float(mm-1)/factory
			    tau_slice(3,i,jj) = (tau(3,i,j)-depth1)/diff  
                
	         endif
	      enddo
		enddo
	    enddo
        
	    nat_slice(i) = jj-1
        
	enddo

    end subroutine
    
!--------------------------------------------------------------------------------------
	subroutine make_mod_tau_unitcell( depth1, depth2, tau, mod_tau, nm, nt, nat, nat2)
	!subroutine make_mod_tau_force() is called if natural depths
      !were NOT chosen.  In this case the current and next depth are
      !required to select which of the tau fit within this slice.
      !Note that a tolerance is employed -- this measures how far
      !beyond the limits an atom can be and still be counted.
      !CAREFUL to realise that false holz can be an issue using this forced
      !slicing.
	use m_precision
    implicit none

	integer(4) nm, nt
	integer(4) nat(nt),nat2(nt)
	real(fp_kind) depth2, depth1, tau(3,nt,nm),mod_tau(3,nt,nm)
	integer(4) i,j,jj
	real(fp_kind) diff, tol
	tol = 1.0e-4_fp_kind

	diff = depth2-depth1


	do i=1,nt
	    jj=1
	    do j=1,nat(i)
	        if( (tau(3,i,j) .lt. (depth2-tol)) .and.(tau(3,i,j) .ge. (depth1-tol)) ) then
	            mod_tau(1,i,jj) = tau(1,i,j)
	            mod_tau(2,i,jj) = tau(2,i,j)
			    mod_tau(3,i,jj) = (tau(3,i,j)-depth1)/diff
			    jj = jj+1
	        elseif(diff.eq.1.0) then
	            mod_tau(1,i,jj) = tau(1,i,j)
	            mod_tau(2,i,jj) = tau(2,i,j)
			    mod_tau(3,i,jj) = (tau(3,i,j)-depth1)/diff  
	        endif
	    enddo
	    nat2(i)=jj-1
	enddo
    
	end subroutine
    
end module
