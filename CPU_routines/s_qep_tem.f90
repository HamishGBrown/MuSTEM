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


subroutine qep_tem
    
    use m_precision, only: fp_kind
	use m_numerical_tools
    use global_variables!, only: nopiy, nopix, ifactory, ifactorx, npixels, nopix_ucell, nopiy_ucell, normalisation, n_cells, ndet, ionization, on_the_fly, ig1, ig2, nt, prop, fz, inverse_sinc, bwl_mat
    use m_lens, only: imaging, imaging_df, pw_illum, probe_df, make_lens_ctf, make_stem_wfn,imaging_ndf,imaging_df,imaging_aberrations,probe_aberrations
	use cufft_wrapper
    use output, only: output_prefix, binary_out, binary_out_unwrap,timing
    use m_slicing
    use m_probe_scan, only: place_probe, probe_initial_position
    use m_tilt, only: tilt_wave_function
    use m_multislice
	use m_string
    use m_potential
    
    implicit none
    
    !dummy variables
    integer(4) :: i_cell, i_slice, i_qep_pass,i,j
    integer(4) :: shifty, shiftx,length,z_indx(1)
    
    !random variables
    integer(4) :: idum
    !real(fp_kind) :: ran1
       
    !probe variables
    complex(fp_kind),dimension(nopiy,nopix) :: psi,psi_initial,psi_out,psi_temp,trans
	complex(fp_kind),dimension(nopiy,nopix,imaging_ndf)::ctf
    complex(fp_kind)::psi_elastic(nopiy,nopix,nz)
	
	
    !output
    real(fp_kind),dimension(nopiy,nopix,nz) :: cbed,total_intensity
	real(fp_kind)::tem_image(nopiy,nopix,nz,imaging_ndf)
    real(fp_kind),dimension(nopiy,nopix) :: image,temp
    
	character*120 ::fnam_df
	integer :: lengthdf
    
    !diagnostic variables
    real(fp_kind) :: intensity, t1, delta
  
    real(fp_kind)::ccd_slice(n_slices)
    real(fp_kind) :: qep_tem_GPU_memory
    
    character*100::filename

    write(*,*) '|----------------------------------|'
	write(*,*) '|      Pre-calculation setup       |'
	write(*,*) '|----------------------------------|'
    write(*,*)

    ccd_slice = relm / (tp * ak * ss_slice(7,:))
    call precalculate_scattering_factors
    !Generally not practical for CPU calculations
    on_the_fly = .false.
    if (on_the_fly) then
        ! call cuda_setup_many_phasegrate
        
    else
        idum = seed_rng()
        
        call make_qep_grates(idum)
        
    endif
    
    call setup_propagators
    
    if(pw_illum) then
    do i=1,imaging_ndf
        call make_lens_ctf(ctf(:,:,i),imaging_df(i),imaging_aberrations)
    enddo
    endif
    
    write(*,*) '|--------------------------------|'
	write(*,*) '|      Calculation running       |'
	write(*,*) '|--------------------------------|'
    write(*,*)
    
	t1 = secnds(0.0)

	if (pw_illum) then
	      psi_initial = 1.0_fp_kind / sqrt(float(npixels))
	else
	      call make_stem_wfn(psi_initial, probe_df(1), probe_initial_position,probe_aberrations)
	endif
    
    call tilt_wave_function(psi_initial)
        
	! Set accumulators to zero
    cbed = 0.0_fp_kind
    psi_elastic = 0.0_fp_kind
    total_intensity = 0.0_fp_kind
    
    do i_qep_pass = 1, n_qep_passes 
    
        ! Reset wavefunction
		psi = psi_initial
        
        do i_cell = 1,maxval(ncells)
	        do i_slice = 1, n_slices
            
                ! Phase grate
				nran = floor(n_qep_grates*ran1(idum)) + 1
				if(on_the_fly) then
					!call make_qep_potential(trans, tau_slice, nat_slice, ss_slice(7,i_slice))
					psi_out = psi*trans
				elseif(quick_shift) then
					shiftx = floor(ifactorx*ran1(idum)) * nopix_ucell
					shifty = floor(ifactory*ran1(idum)) * nopiy_ucell
					trans = cshift(cshift(qep_grates(:,:,nran,i_slice),shifty,dim=1),shiftx,dim=2)
					psi_out = psi*trans
				elseif(phase_ramp_shift) then                       !randomly shift phase grate
					shiftx = floor(ifactorx*ran1(idum)) + 1
					shifty = floor(ifactory*ran1(idum)) + 1
					call phase_shift_array(qep_grates(:,:,nran,i_slice),trans,shift_arrayy,shift_arrayx)
					psi_out = psi*trans
				else
					psi_out = psi*qep_grates(:,:,nran,i_slice)
				endif
                
				!propagate
                call fft2(nopiy,nopix,psi_out,nopiy,psi,nopiy)
				psi_out = prop(:,:,i_slice)*psi
				call ifft2(nopiy,nopix,psi_out,nopiy,psi,nopiy)


            enddo ! End loop over slices
			
				!If this thickness is an output thickness then accumulate relevent TEM images and diffraction patterns
				if (any(i_cell==ncells)) then
					
					!Transform into diffraction space
					call fft2(nopiy,nopix,psi,nopiy,psi_out,nopiy)
					temp = abs(psi_out)**2
					z_indx = minloc(abs(ncells-i_cell))

					! Accumulate elastic wave function
					psi_elastic(:,:,z_indx(1)) = psi_elastic(:,:,z_indx(1)) + psi

					! Accumulate exit surface intensity
					total_intensity(:,:,z_indx(1)) = total_intensity(:,:,z_indx(1)) + abs(psi)**2
					
					! Accumulate diffaction pattern
					temp=abs(psi_out)**2
					cbed(:,:,z_indx(1)) = cbed(:,:,z_indx(1)) + temp
					
					! Accumulate image
                    if(pw_illum) then
					do i=1,imaging_ndf
						psi_temp = ctf(:,:,i)*psi_out
						call ifft2(nopiy,nopix,psi_temp,nopiy,psi_temp,nopiy)
						tem_image(:,:,z_indx(1),i) = tem_image(:,:,z_indx(1),i)+abs(psi_temp)**2
                    enddo
                    endif

				endif
        enddo ! End loop over cells
        
        intensity = sum(abs(psi)**2)
		
#ifdef gpu       
    write(6,900,advance='no') achar(13), i_qep_pass, n_qep_passes, intensity
    900 format(a1, 1x, 'QEP pass:', i4, '/', i4, ' Intensity: ', f8.3)	
#else
        write(6,900) i_qep_pass, n_qep_passes, intensity
900     format(1h+,  'QEP pass:', i4, '/', i4, ' Intensity: ', f8.3)	
#endif
	enddo ! End loop over QEP passes
   
    delta = secnds(t1)
    
    write(*,*)
    write(*,*)
    write(*,*) 'Calculation is finished.'
    write(*,*)
    write(*,*) 'Time elapsed: ', delta, ' seconds.'
    write(*,*)  
    
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
    
    cbed = cbed / n_qep_passes
    total_intensity = total_intensity / n_qep_passes
    psi_elastic = psi_elastic / n_qep_passes
    tem_image = tem_image / n_qep_passes
    
    length = ceiling(log10(maxval(zarray)))
	lengthdf = ceiling(log10(maxval(abs(imaging_df))))
	if(any(imaging_df<0)) lengthdf = lengthdf+1
   
	do i=1,nz
		filename = trim(adjustl(output_prefix))
		if(nz>1) filename=trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A'
		
        
        if(.not.output_thermal) then
		call binary_out_unwrap(nopiy, nopix, cbed(:,:,i), trim(filename)//'_DiffPlane')
        call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_Intensity')
        else
        call binary_out_unwrap(nopiy, nopix, cbed(:,:,i), trim(filename)//'_DiffPlaneTotal')
		
		call fft2(nopiy, nopix, psi_elastic(:,:,i), nopiy, psi, nopiy)
		image = abs(psi)**2
		call binary_out_unwrap(nopiy, nopix, image, trim(filename)//'_DiffPlaneElastic')
		
		image = cbed(:,:,i) - image
		call binary_out_unwrap(nopiy, nopix, image, trim(filename)//'_DiffPlaneTDS')
		
		call binary_out(nopiy, nopix, abs(psi_elastic(:,:,i))**2, trim(filename)//'_ExitSurface_IntensityElastic')
		call binary_out(nopiy, nopix, atan2(imag(psi_elastic(:,:,i)), real(psi_elastic(:,:,i))), trim(filename)//'_ExitSurface_PhaseElastic')
		call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_IntensityTotal')
		
		total_intensity(:,:,i) = total_intensity(:,:,i) - abs(psi_elastic(:,:,i))**2
		call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_IntensityTDS')
        endif
            if (pw_illum) then
			do j=1,imaging_ndf
				! Elastic image
				call fft2 (nopiy, nopix, psi_elastic(:,:,i), nopiy, psi, nopiy)
				psi = psi * ctf(:,:,j)
				call ifft2 (nopiy, nopix, psi, nopiy, psi, nopiy)
				image = abs(psi)**2
				
				if(imaging_ndf>1) fnam_df = trim(adjustl(filename))//'_Defocus_'//zero_padded_int(int(imaging_df(j)),lengthdf)//'_Ang'
                if(output_thermal) then
				call binary_out(nopiy, nopix, image, trim(fnam_df)//'_Image_Elastic')
				  
				! Inelastic image
				image = tem_image(:,:,i,j) - image
				call binary_out(nopiy, nopix, image, trim(fnam_df)//'_Image_TDS')
                
				! Total image
				call binary_out(nopiy, nopix, tem_image(:,:,i,j), trim(fnam_df)//'_Image_Total')
                else
                call binary_out(nopiy, nopix, tem_image(:,:,i,j), trim(fnam_df)//'_Image')
                endif
            enddo
            endif
	
	enddo
end subroutine qep_tem