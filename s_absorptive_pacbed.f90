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

    subroutine absorptive_pacbed
      
        use global_variables, only: nopiy, nopix, n_cells, prop
        use m_lens, only: make_stem_wfn,probe_df
        use m_absorption, only: calculate_absorption_mu, transf_absorptive
        use cufft_wrapper, only: fft2, ifft2
        use m_precision, only: fp_kind
        use output, only: output_prefix, binary_out_unwrap
        use m_slicing, only: n_slices
        use m_probe_scan, only: nysample, nxsample, probe_positions
        use m_tilt, only: tilt_wave_function
        use m_multislice, only: make_absorptive_grates, setup_propagators, setup_propagators
        use m_potential, only: precalculate_scattering_factors
		use m_string
		use m_user_input
        
        implicit none
      
        !dummy variables
        integer(4) ::  i_cell, i_slice, ny, nx, idum

        !probe variables
        complex(fp_kind) :: psi(nopiy,nopix)
      
        !output/detectors
        real(fp_kind),dimension(nopiy, nopix) :: pacbed_pattern,fourDSTEM_pattern

        !diagnostic variables
        real(fp_kind) :: intensity, t1, delta
      
        !output variables
        character(120) :: filename
		logical::fourdstem

        write(*,*) '|----------------------------------|'
	    write(*,*) '|      Pre-calculation setup       |'
	    write(*,*) '|----------------------------------|'
        write(*,*)


        call calculate_absorption_mu        

        call precalculate_scattering_factors()
    
        call make_absorptive_grates
        
        call setup_propagators()

		
        write(*,*) '|--------------------------------|'
	    write(*,*) '|      Calculation running       |'
	    write(*,*) '|--------------------------------|'
        write(*,*)
      
        t1 = secnds(0.0)
      
        pacbed_pattern = 0.0_fp_kind

        intensity = 1.0_fp_kind
        write(*,*) 'Output diffraction patterns for each scan position?'
		call get_input('<1> Diffraction pattern for each probe position',idum)
		fourDSTEM = idum == 1
        do ny = 1, nysample
        do nx = 1, nxsample
            write(6, 901, advance='no') achar(13), ny, nysample, nx, nxsample, intensity
901         format(a1, 1x, 'y:', i3, '/', i3, ' x:', i3, '/', i3, ' Intensity: ', f12.6)
      
            call make_stem_wfn(psi, probe_df, probe_positions(:,ny,nx))
            
            call tilt_wave_function(psi)
            
            do i_cell = 1, n_cells
                  do i_slice = 1, n_slices
	                    psi = psi*transf_absorptive(:,:,i_slice)
                        
                        call fft2(nopiy, nopix, psi, nopiy, psi, nopiy)
                        psi = prop(:,:,i_slice) * psi
                        call ifft2(nopiy, nopix, psi, nopiy, psi, nopiy)
                  enddo
            enddo
            
            call fft2(nopiy, nopix, psi, nopiy, psi, nopiy)
			fourDSTEM_pattern = abs(psi)**2
            pacbed_pattern = pacbed_pattern + fourDSTEM_pattern

			!Output 4D STEM diffraction pattern
			if(fourDSTEM) then
					filename = trim(adjustl(output_prefix)) //'_pp_'//to_string(nx)//'_'//to_string(ny)//'_abs_Diffraction_pattern'
					call binary_out_unwrap(nopiy, nopix, fourDSTEM_pattern, filename,write_to_screen=.false.)
			endif
            
            intensity = sum(fourDSTEM_pattern)
        enddo
        enddo
      
        delta = secnds(t1)
      
        write(*,*)                 
        write(*,*)
        
        write(*,*) 'Calculation is finished.'
        write(*,*) 
        write(*,*) 'Time elapsed ', delta, ' seconds.'
        
        open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
        write(9834, '(a, g, a, /)') 'The multislice calculation took ', delta, 'seconds.'
        close(9834)
    
      
        if (fp_kind.eq.8) then
            write(*,*) 'The following files were outputted (as 64-bit big-endian floating point):'
	    else
            write(*,*) 'The following files were outputted (as 32-bit big-endian floating point):'
	    endif
        write(*,*)
    
        filename = trim(adjustl(output_prefix)) // '_PACBED_Pattern'
        call binary_out_unwrap(nopiy,nopix,pacbed_pattern,filename)
      
        deallocate(transf_absorptive,prop) !deallocate large arrays

    end subroutine absorptive_pacbed
      