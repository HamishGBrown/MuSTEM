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
      
        use global_variables, only: nopiy, nopix, n_cells, prop,zarray,nz,ncells
        use m_lens, only: make_stem_wfn,probe_df,probe_aberrations
        use m_absorption, only: calculate_absorption_mu, transf_absorptive
        use cufft_wrapper, only: fft2, ifft2
        use m_precision, only: fp_kind
        use output, only: output_prefix, binary_out_unwrap,timing
        use m_slicing, only: n_slices
        use m_probe_scan, only: nysample, nxsample, probe_positions
        use m_tilt, only: tilt_wave_function
        use m_multislice, only: make_absorptive_grates, setup_propagators, setup_propagators
        use m_potential, only: precalculate_scattering_factors
		use m_string
		use m_user_input
        
        implicit none
      
        !dummy variables
        integer(4) ::  i_cell, i_slice, ny, nx, idum,i,z_indx(1),length,nopiyout,nopixout

        !probe variables
        complex(fp_kind),dimension(nopiy,nopix) :: psi,psi_
      
        !output/detectors
        real(fp_kind) :: pacbed_pattern(nopiy, nopix,nz),fourDSTEM_pattern(nopiy, nopix)

        !diagnostic variables
        real(fp_kind) :: intensity, t1, delta
      
        !output variables
        character(120) :: filename
		logical::fourdstem

        write(*,*) '|----------------------------------|'
	    write(*,*) '|      Pre-calculation setup       |'
	    write(*,*) '|----------------------------------|'
        write(*,*)    

        call precalculate_scattering_factors()
    
        call make_absorptive_grates
        
        call setup_propagators()
        call fourD_STEM_options(fourdSTEM,nopiyout,nopixout,nopiy,nopix)
		
        write(*,*) '|--------------------------------|'
	    write(*,*) '|      Calculation running       |'
	    write(*,*) '|--------------------------------|'
        write(*,*)
      
        t1 = secnds(0.0)
        
        pacbed_pattern = 0.0_fp_kind

        intensity = 1.0_fp_kind
        
        do ny = 1, nysample
        do nx = 1, nxsample
#ifdef GPU                    
            write(6, 901, advance='no') achar(13), ny, nysample, nx, nxsample, intensity
901         format(a1, 1x, 'y:', i3, '/', i3, ' x:', i3, '/', i3, ' Intensity: ', f12.6)
#else
        write(6,900) ny, nysample, nx, nxsample, intensity
900     format(1h+, 1x, 'y:', i3, '/', i3, ' x:', i3, '/', i3, ' Intensity: ', f12.6)
#endif      
            call make_stem_wfn(psi, probe_df(1), probe_positions(:,ny,nx),probe_aberrations)
            
            call tilt_wave_function(psi)
            
			do i_cell = 1,maxval(ncells)
                  do i_slice = 1, n_slices
	                    psi = psi*transf_absorptive(:,:,i_slice)
                        
                        call fft2(nopiy, nopix, psi, nopiy, psi, nopiy)
                        psi = prop(:,:,i_slice) * psi
                        call ifft2(nopiy, nopix, psi, nopiy, psi, nopiy)
                  enddo

				!If this thickness corresponds to any of the output values then output images
				if (any(i_cell==ncells)) then
					call fft2(nopiy, nopix, psi, nopiy, psi_, nopiy)
					fourDSTEM_pattern = abs(psi_)**2
					z_indx = minloc(abs(ncells-i_cell))
					pacbed_pattern(:,:,z_indx(1)) = pacbed_pattern(:,:,z_indx(1)) + fourDSTEM_pattern

					!Output 4D STEM diffraction pattern
					if(fourDSTEM) then
							filename = trim(adjustl(output_prefix))
							if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(z_indx(1))))//'_A'
							filename = trim(adjustl(filename))//'_pp_'//to_string(nx)//'_'//to_string(ny)//'_abs_Diffraction_pattern'
							call binary_out_unwrap(nopiy, nopix, fourDSTEM_pattern, filename,write_to_screen=.false.,nopiyout=nopiyout,nopixout=nopixout)
					endif			
				endif
            enddo
            
            intensity = sum(fourDSTEM_pattern)
        enddo
        enddo
      
        delta = secnds(t1)
      
        write(*,*)                 
        write(*,*)
        
        write(*,*) 'Calculation is finished.'
        write(*,*) 
        write(*,*) 'Time elapsed ', delta, ' seconds.'
        
		if(timing) then
			open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
			write(9834, '(a, g, a, /)') 'The multislice calculation took ', delta, 'seconds.'
			close(9834)
		endif
    
      
        if (fp_kind.eq.8) then
            write(*,*) 'The following files were outputted (as 64-bit big-endian floating point):'
	    else
            write(*,*) 'The following files were outputted (as 32-bit big-endian floating point):'
	    endif
        write(*,*)
		
		length = ceiling(log10(maxval(zarray)))
		do i=1,nz
			filename = trim(adjustl(output_prefix))
			if(nz>1) filename = trim(adjustl(output_prefix))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A'
			filename = trim(adjustl(filename))//'_PACBED_Pattern'
			call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i),filename)
		enddo
    
        ! filename = trim(adjustl(output_prefix)) // '_PACBED_Pattern'
        ! call binary_out_unwrap(nopiy,nopix,pacbed_pattern,filename)
      
        deallocate(transf_absorptive,prop) !deallocate large arrays

    end subroutine absorptive_pacbed
      