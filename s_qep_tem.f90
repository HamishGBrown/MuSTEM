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

function qep_tem_GPU_memory() result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, on_the_fly
    use m_lens, only: imaging
    use m_slicing, only: n_slices
    use m_qep, only: n_qep_grates, quick_shift, phase_ramp_shift
    
    implicit none
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 2*(3 + n_slices + n_slices*n_qep_grates) + 3
    
    if (imaging) array_count = array_count + 2 + 1
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) array_count = array_count + 2
    if (phase_ramp_shift) array_count = array_count + 2
    
    required_memory = array_count * array_size
    
    if (phase_ramp_shift) required_memory = required_memory + 8.0_fp_kind*(nopiy*ifactory + nopix*ifactorx)
    
end function
    


subroutine qep_tem
    
    use m_precision, only: fp_kind
	use m_numerical_tools
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, npixels, nopix_ucell, nopiy_ucell, normalisation, n_cells, ndet, ionization, on_the_fly, ig1, ig2, nt, prop, fz, inverse_sinc, bwl_mat
    use m_lens, only: imaging, imaging_df, pw_illum, probe_df, make_ctf, make_stem_wfn
    use m_qep, only: n_qep_passes, n_qep_grates, seed_rng, qep_grates, quick_shift, phase_ramp_shift, shift_arrayx, shift_arrayy, nran
	use cufft_wrapper, only: cufft_z2z, cufft_c2c, cufft_inverse, cufft_forward, cufftPlan, cufftExec, fft2, ifft2
    use cuda_array_library, only: cuda_multiplication, cuda_addition, cuda_mod, blocks, threads
    use cudafor, only: dim3
    use cuda_ms, only: get_sum
    use output, only: output_prefix, binary_out, binary_out_unwrap
    use cuda_potential, only: cuda_setup_many_phasegrate, ccd_slice_array, cuda_fph_make_potential
    use m_slicing, only: n_slices, nat_slice, tau_slice, prop_distance
    use cuda_setup, only: GPU_memory_message
    use m_probe_scan, only: place_probe, probe_initial_position
    use m_tilt, only: tilt_wave_function
    use m_multislice, only: make_qep_grates, setup_propagators
    use m_potential, only: precalculate_scattering_factors
    
    implicit none
    
    !dummy variables
    integer(4) :: i_cell, i_slice, i_qep_pass
    integer(4) :: shifty, shiftx
    
    !random variables
    integer(4) :: idum
    !real(fp_kind) :: ran1
       
    !probe variables
    complex(fp_kind) :: psi(nopiy,nopix)
    complex(fp_kind) :: psi_initial(nopiy,nopix)
    complex(fp_kind) :: psi_elastic(nopiy,nopix)
    complex(fp_kind) :: ctf(nopiy,nopix)
    
    !output
    real(fp_kind) :: cbed(nopiy,nopix)
    real(fp_kind) :: image(nopiy,nopix)
    real(fp_kind) :: tem_image(nopiy,nopix)
    real(fp_kind) :: total_intensity(nopiy,nopix)
    
    !diagnostic variables
    real(fp_kind) :: intensity, t1, delta
    
    !device variables
	integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d, psi_elastic_d, psi_initial_d
	complex(fp_kind),device,allocatable :: prop_d(:,:,:), transf_d(:,:,:,:), ctf_d(:,:), shift_arrayx_d(:,:), shift_arrayy_d(:,:), shift_array_d(:,:), trans_d(:,:)
    real(fp_kind),device,allocatable :: tem_image_d(:,:)
    real(fp_kind),device,dimension(nopiy,nopix) :: cbed_d, temp_image_d, total_intensity_d
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d, inverse_sinc_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d

    real(fp_kind) :: qep_tem_GPU_memory
    
    
    
    call GPU_memory_message(qep_tem_GPU_memory(), on_the_fly)
    
    
    
    write(*,*) '|----------------------------------|'
	write(*,*) '|      Pre-calculation setup       |'
	write(*,*) '|----------------------------------|'
    write(*,*)

    call precalculate_scattering_factors

    if (on_the_fly) then
        call cuda_setup_many_phasegrate
        
    else
        idum = seed_rng()
        
        call make_qep_grates(idum)
        
    endif
    
    call setup_propagators
    if (imaging) then
          call make_ctf(ctf,imaging_df)
          allocate(ctf_d(nopiy,nopix))
          ctf_d = ctf
          
          allocate(tem_image_d(nopiy,nopix))
          tem_image_d = 0.0_fp_kind
    endif

    
    
    write(*,*) '|--------------------------------|'
	write(*,*) '|      Calculation running       |'
	write(*,*) '|--------------------------------|'
    write(*,*)

    if (fp_kind.eq.8)then
	    call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
	else
        call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
    
	allocate(prop_d(nopiy,nopix,n_slices))
	prop_d = prop

    if (on_the_fly) then
        if(allocated(bwl_mat_d)) deallocate(bwl_mat_d)
        if(allocated(inverse_sinc_d)) deallocate(inverse_sinc_d)
        if(allocated(fz_d)) deallocate(fz_d)
        allocate(bwl_mat_d(nopiy,nopix))
        allocate(inverse_sinc_d(nopiy,nopix))
        allocate(fz_d(nopiy,nopix,nt))
        fz_d = fz
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
    else
        allocate(transf_d(nopiy,nopix,n_qep_grates,n_slices))
        transf_d = qep_grates
        
        if (phase_ramp_shift) then
            allocate(shift_array_d(nopiy,nopix))
            allocate(shift_arrayx_d(nopix,ifactorx))
            allocate(shift_arrayy_d(nopiy,ifactory))
            shift_arrayx_d=shift_arrayx
            shift_arrayy_d=shift_arrayy
        endif
    endif
    
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) then
        allocate(trans_d(nopiy,nopix))
    endif
    
	t1 = secnds(0.0)

	if (pw_illum) then
	      psi_initial = 1.0_fp_kind / sqrt(float(npixels))
	else
	      call make_stem_wfn(psi_initial, probe_df, probe_initial_position)
	endif
    
    call tilt_wave_function(psi_initial)
    
    psi_initial_d = psi_initial
        
	! Set accumulators to zero
    cbed_d = 0.0_fp_kind
    psi_elastic_d = 0.0_fp_kind
    total_intensity_d = 0.0_fp_kind
    
    do i_qep_pass = 1, n_qep_passes 
    
        ! Reset wavefunction
        psi_d = psi_initial_d
        
        do i_cell = 1, n_cells
	        do i_slice = 1, n_slices
            
                ! Transmit
                if (on_the_fly) then
                    call cuda_fph_make_potential(trans_d,ccd_slice_array(i_slice),tau_slice,nat_slice(:,i_slice),i_slice,prop_distance(i_slice),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
                    call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_d,1.0_fp_kind,nopiy,nopix)
                elseif (quick_shift) then
				    nran = floor(n_qep_grates*ran1(idum)) + 1
                    shiftx = floor(ifactorx*ran1(idum)) * nopix_ucell
                    shifty = floor(ifactory*ran1(idum)) * nopiy_ucell
                    call cuda_cshift<<<blocks,threads>>>(transf_d(:,:,nran,i_slice),trans_d,nopiy,nopix,shifty,shiftx)
                    call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_d,1.0_fp_kind,nopiy,nopix)
                elseif (phase_ramp_shift) then
				    nran = floor(n_qep_grates*ran1(idum)) + 1
                    shiftx = floor(ifactorx*ran1(idum)) + 1
                    shifty = floor(ifactory*ran1(idum)) + 1
                    call cuda_make_shift_array<<<blocks,threads>>>(shift_array_d,shift_arrayy_d(:,shifty),shift_arrayx_d(:,shiftx),nopiy,nopix)
                    call cuda_multiplication<<<blocks,threads>>>(transf_d(:,:,nran,i_slice),shift_array_d, trans_d,1.0_fp_kind,nopiy,nopix)
                    call cufftExec(plan,trans_d,trans_d,CUFFT_INVERSE)
                    call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_d,sqrt(normalisation),nopiy,nopix)
                else
				    nran = floor(n_qep_grates*ran1(idum)) + 1
                    call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d(:,:,nran,i_slice), psi_d,1.0_fp_kind,nopiy,nopix)
                endif
                
                ! Propagate
				call cufftExec(plan, psi_d, psi_d, CUFFT_FORWARD)
                call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d(:,:,i_slice), psi_d, normalisation, nopiy, nopix)
                call cufftExec(plan, psi_d, psi_d, CUFFT_INVERSE)


            enddo ! End loop over slices
        enddo ! End loop over cells
        
        intensity = get_sum(psi_d)
		write(6,900,advance='no') achar(13), i_qep_pass, n_qep_passes, intensity
    900 format(a1, 1x, 'QEP pass:', i4, '/', i4, ' Intensity: ', f8.3)	
    
        ! Accumulate elastic wave function
        call cuda_addition<<<blocks,threads>>>(psi_elastic_d, psi_d, psi_elastic_d, 1.0_fp_kind, nopiy, nopix)
        
        ! Accumulate exit surface intensity
        call cuda_mod<<<blocks,threads>>>(psi_d, temp_image_d, 1.0_fp_kind, nopiy, nopix)
        call cuda_addition<<<blocks,threads>>>(total_intensity_d, temp_image_d, total_intensity_d, 1.0_fp_kind, nopiy, nopix)
        
        ! Accumulate diffraction pattern
        call cufftExec(plan, psi_d, psi_d, CUFFT_FORWARD)
        call cuda_mod<<<blocks,threads>>>(psi_d, temp_image_d, normalisation, nopiy, nopix)
        call cuda_addition<<<blocks,threads>>>(cbed_d, temp_image_d, cbed_d, 1.0_fp_kind, nopiy, nopix)
        
        ! Accumulate image
        if (imaging) then
              call cuda_multiplication<<<blocks,threads>>>(psi_d, ctf_d, psi_d, sqrt(normalisation), nopiy, nopix)
              call cufftExec(plan, psi_d, psi_d, CUFFT_INVERSE)
              call cuda_mod<<<blocks,threads>>>(psi_d, temp_image_d, normalisation, nopiy, nopix)
              call cuda_addition<<<blocks,threads>>>(tem_image_d, temp_image_d, tem_image_d, 1.0_fp_kind, nopiy, nopix)
        endif

	enddo ! End loop over QEP passes
    
    cbed = cbed_d
    psi_elastic = psi_elastic_d
    total_intensity = total_intensity_d
    if (imaging) then
        tem_image = tem_image_d
    endif      
    
    delta = secnds(t1)
    
    write(*,*)
    write(*,*)
    write(*,*) 'Calculation is finished.'
    write(*,*)
    write(*,*) 'Time elapsed: ', delta, ' seconds.'
    write(*,*)  
    
    open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
    write(9834, '(a, g, a, /)') 'The multislice calculation took ', delta, 'seconds.'
    close(9834)
    
    if (fp_kind.eq.8) then
        write(*,*) 'The following files were outputted (as 64-bit big-endian floating point):'
	else
        write(*,*) 'The following files were outputted (as 32-bit big-endian floating point):'
	endif
    write(*,*)
    
    cbed = cbed / n_qep_passes
    total_intensity = total_intensity / n_qep_passes
    psi_elastic = psi_elastic / n_qep_passes
    if (imaging) then
        tem_image = tem_image / n_qep_passes
    endif      
    
    
    call binary_out_unwrap(nopiy, nopix, cbed, trim(output_prefix)//'_DiffPlaneTotal')
    
    call fft2(nopiy, nopix, psi_elastic, nopiy, psi, nopiy)
    image = abs(psi)**2
    call binary_out_unwrap(nopiy, nopix, image, trim(output_prefix)//'_DiffPlaneElastic')
    
    image = cbed - image
    call binary_out_unwrap(nopiy, nopix, image, trim(output_prefix)//'_DiffPlaneTDS')
    
    call binary_out(nopiy, nopix, abs(psi_elastic)**2, trim(output_prefix)//'_ExitSurface_IntensityElastic')
    call binary_out(nopiy, nopix, atan2(imag(psi_elastic), real(psi_elastic)), trim(output_prefix)//'_ExitSurface_PhaseElastic')
    call binary_out(nopiy, nopix, total_intensity, trim(output_prefix)//'_ExitSurface_IntensityTotal')
    
    total_intensity = total_intensity - abs(psi_elastic)**2
    call binary_out(nopiy, nopix, total_intensity, trim(output_prefix)//'_ExitSurface_IntensityTDS')

    if (imaging) then          
        ! Elastic image
        call fft2 (nopiy, nopix, psi_elastic, nopiy, psi_elastic, nopiy)
        psi_elastic = psi_elastic * ctf
        call ifft2 (nopiy, nopix, psi_elastic, nopiy, psi_elastic, nopiy)
        image = abs(psi_elastic)**2
        call binary_out(nopiy, nopix, image, trim(output_prefix)//'_Image_Elastic')
          
        ! Inelastic image
        image = tem_image - image
        call binary_out(nopiy, nopix, image, trim(output_prefix)//'_Image_TDS')
          
        ! Total image
        call binary_out(nopiy, nopix, tem_image, trim(output_prefix)//'_Image_Total')
    endif

end subroutine qep_tem
