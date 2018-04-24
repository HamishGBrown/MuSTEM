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

function qep_stem_GPU_memory() result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, on_the_fly, ndet, ionization, eels
    use m_lens, only: imaging
    use m_slicing, only: n_slices
    use m_qep, only: n_qep_grates, quick_shift, phase_ramp_shift
    
    implicit none
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 2*(4 + n_slices + n_slices*n_qep_grates) + 2 + ndet
        
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) array_count = array_count + 2
    if (phase_ramp_shift) array_count = array_count + 2
    if (ionization) array_count = array_count + 1 + n_slices
    if (eels) array_count = array_count + 1
    
    ! Temporary array used in cuda_stem_detector()
    array_count = array_count + 1
    
    required_memory = array_count * array_size
    
    if (phase_ramp_shift) required_memory = required_memory + 8.0_fp_kind*(nopiy*ifactory + nopix*ifactorx)
    
end function

    

subroutine qep_stem
    
    use global_variables
    use m_lens
    use m_user_input
    use m_qep
    use m_precision
    use output
	use cufft_wrapper
    use local_ionization
    use m_slicing
    use m_probe_scan, only: nysample, nxsample, probe_positions, scan_quarter
    use m_tilt, only: tilt_wave_function
    use m_string, only: to_string
    use m_multislice!, only: make_qep_grates, output_probe_intensity, output_cell_list, cell_map, output_thickness_list, setup_propagators
    use m_potential
    
    implicit none
    
    !dummy variables
    integer(4) ::  i,j,l,m,i_qep_pass,iz
    integer(4) ::  count
    integer(4) ::  shifty,shiftx
    integer(4) ::  ny,nx,i_df,idet
    integer(4) ::  total_slices
    
    !random variables
    integer(4) :: idum,nsliceoutput,z_indx(1),k,ii,jj
    !real(fp_kind) :: ran1
    complex(fp_kind) :: temp_transf(nopiy,nopix)
       
    !probe variables
    complex(fp_kind),dimension(nopiy,nopix) :: psi,trans,psi_initial,psi_out
    complex(fp_kind),dimension(nopiy,nopix,nz)::psi_elastic
     
    !output
    real(fp_kind),dimension(nopiy,nopix) :: image,psi_intensity,inelastic_potential,temp,temp_image
    real(fp_kind),dimension(nopiy,nopix,nz) :: cbed
    real(fp_kind),allocatable::ion_image(:,:,:)
	
    !STEM image variables
    real(fp_kind) :: masks(nopiy,nopix,ndet)                       !detector masks
    real(fp_kind),dimension(nysample,nxsample,probe_ndf,ndet,nz) :: stem_image,stem_elastic_image
    real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz) :: stem_inelastic_image,eels_correction_image
    real(fp_kind),allocatable::stem_ion_image(:,:,:,:,:)
    !diagnostic variables
    real(fp_kind) :: intensity
    real(fp_kind) :: t1, delta
    
    !output variables
    character(120) :: fnam, fnam_temp, fnam_det

	!Probe depth output
	!real(fp_kind),allocatable :: probe_intensity_output(:,:,:)

    real(fp_kind),allocatable :: probe_intensity(:,:,:)
    character(1024) :: filename
    
    real(fp_kind)::qep_stem_GPU_memory
    
    write(*,*) '|----------------------------------|'
	write(*,*) '|      Pre-calculation setup       |'
	write(*,*) '|----------------------------------|'
    write(*,*)
            
    !Generally not practical for CPU calculations
    on_the_fly = .false.
    idum = seed_rng()
    if (.not.on_the_fly) then
        call make_qep_grates(idum)
        call make_local_inelastic_potentials()          !setup the REAL space inelastic potential (ionization and adf) for QUEP ADF is disabled
    endif
    
 	call setup_propagators
    
    do i = 1, ndet
        call make_detector_mask(inner(i),outer(i),masks(:,:,i))
    enddo

    
	t1 = secnds(0.0)
    
    
    write(*,*) '|--------------------------------|'
	write(*,*) '|      Calculation running       |'
	write(*,*) '|--------------------------------|'
    write(*,*)
       
    if (output_probe_intensity) allocate(probe_intensity(nopiy,nopix,size(output_thickness_list)))
    
    stem_image = 0.0_fp_kind
    stem_elastic_image = 0.0_fp_kind
    if (ionization) then 
        allocate(ion_image(nopiy,nopix,num_ionizations))
        allocate(stem_ion_image(nysample,nxsample,probe_ndf,nz,num_ionizations))
        stem_ion_image = 0.0_fp_kind
		if (.not.EDX) eels_correction_image = 0.0_fp_kind
    endif
    
    
    intensity = 1.0d0
    
	do i_df = 1, probe_ndf
	do ny = 1, nysample
    do nx = 1, nxsample
        
#ifdef gpu       
    write(6,903,ADVANCE='NO') achar(13), i_df, probe_ndf, ny, nysample, nx, nxsample, intensity
903     format(a1,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3, ' (to monitor BWL)')	
#else
        write(6,900) i_df, probe_ndf, ny, nysample, nx, nxsample, intensity
900     format(1h+,1x,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3, ' (to monitor BWL)')	
#endif        
    
        if (scan_quarter) then
            if (probe_positions(1,ny,nx).gt.0.501_fp_kind .or. probe_positions(2,ny,nx).gt.0.501_fp_kind) cycle
        endif
    

        cbed=0.0_fp_kind
        
        psi_elastic=0.0_fp_kind
        call make_stem_wfn(psi_initial,probe_df(i_df),probe_positions(:,ny,nx),probe_aberrations)
        
        call tilt_wave_function(psi_initial)
        
        if (output_probe_intensity) probe_intensity = 0.0_fp_kind
        
        do i_qep_pass = 1, n_qep_passes 
			if (ionization) ion_image=0.0_fp_kind
            ! Reset wavefunction
            psi = psi_initial
            
            do i = 1,maxval(ncells)
              
	            do j = 1, n_slices
                    ! Accumulate ionization cross section
					if(ionization) then
						psi_intensity = abs(psi)**2
                        do ii=1,num_ionizations
                            if(on_the_fly) then                        
                                inelastic_potential= make_ion_potential(ionization_mu(:,:,ii),tau_slice(:,atm_indices(ii),:,j),nat_slice(atm_indices(ii),j),ss_slice(7,j))
						        temp = psi_intensity *inelastic_potential*prop_distance(j)
                            else
						        temp = psi_intensity * ionization_potential(:,:,ii,j) * prop_distance(j)
                            endif
					        ion_image(:,:,ii) = temp+ion_image(:,:,ii)
                        enddo
					endif
                    
                   ! Phase grate
				    nran = floor(n_qep_grates*ran1(idum)) + 1
                    if(on_the_fly) then
						call make_qep_potential(trans, tau_slice, nat_slice, ss_slice(7,j))
						psi = psi*trans
                    elseif(quick_shift) then
                        shiftx = floor(ifactorx*ran1(idum)) * nopix_ucell
                        shifty = floor(ifactory*ran1(idum)) * nopiy_ucell
						trans = cshift(cshift(qep_grates(:,:,nran,j),shifty,dim=1),shiftx,dim=2)
						psi = psi*trans
                    elseif(phase_ramp_shift) then                       !randomly shift phase grate
                        shiftx = floor(ifactorx*ran1(idum)) + 1
                        shifty = floor(ifactory*ran1(idum)) + 1
						call phase_shift_array(qep_grates(:,:,nran,j),trans,shift_arrayy,shift_arrayx)
                        psi = psi*trans
                    else
						psi = psi*qep_grates(:,:,nran,j)
                    endif
					

					call fft2(nopiy,nopix,psi,nopiy,psi,nopiy)
					psi = psi*prop(:,:,j)
					call ifft2(nopiy,nopix,psi,nopiy,psi,nopiy)
                    
					if (output_probe_intensity) then
						k = (i-1)*n_slices+j
						if (output_cell_list(k)) then
							probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
						endif
					endif
		        enddo ! End loop over slices
                    
				
                !If this thickness corresponds to any of the output values then accumulate diffration pattern
				if (any(i==ncells)) then
					z_indx = minloc(abs(ncells-i))

					! Accumulate elastic wave function
					psi_elastic(:,:,z_indx(1)) = psi_elastic(:,:,z_indx(1)) + psi

					!Transform into diffraction space
					call fft2(nopiy,nopix,psi,nopiy,psi_out,nopiy)
					! Accumulate diffaction pattern
					temp = abs(psi_out)**2
					cbed(:,:,z_indx(1)) = cbed(:,:,z_indx(1)) + temp

					if(ionization) then
                        stem_ion_image(ny,nx,i_df,z_indx(1),:) = stem_ion_image(ny,nx,i_df,z_indx(1),:)+ sum(sum(ion_image,dim=2),dim=1)
                        if(.not.EDX) eels_correction_image(ny,nx,i_df,z_indx(1)) = sum(temp*eels_correction_detector)
                    endif
				endif
        

                
            enddo ! End loop over cells
            
	    enddo ! End loop over QEP passes
        !call binary_out(nopiy,nopix,abs(psi)**2,'psi_final')
        intensity = sum(abs(psi)**2)
        
		do iz=1,nz
			! Integrate the diffraction pattern
			do idet = 1, ndet
				stem_image(ny,nx,i_df,idet,iz) = sum(cbed(:,:,iz)*masks(:,:,idet))
			enddo
        
        
			if (eels) eels_correction_image(ny,nx,i_df,iz) = sum(abs(cbed(:,:,iz))**2*masks(:,:,idet))
        
			! Integrate the elastic diffraction pattern
			call fft2(nopiy,nopix,psi_elastic(:,:,iz),nopiy,psi_out,nopiy)
			do idet = 1, ndet
				stem_elastic_image(ny,nx,i_df,idet,iz)= sum(abs(psi_out)**2*masks(:,:,idet))
			enddo
		enddo

        if (output_probe_intensity) call probe_intensity_to_file(probe_intensity,i_df,ny,nx,n_qep_passes)
	enddo ! End loop over x probe positions
	enddo ! End loop over y probe positions
	enddo ! End loop over defocus series
	
    ! QEP normalisation
    
    if(ionization) stem_ion_image = stem_ion_image/float(n_qep_passes)
    
    if(EELS) eels_correction_image = eels_correction_image/float(n_qep_passes)

    if (ndet.gt.0) then
	    stem_image = stem_image/float(n_qep_passes)
        stem_elastic_image = stem_elastic_image/float(n_qep_passes*n_qep_passes)
    endif

    delta = secnds(t1)
    
    write(*,*) 
    write(*,*) 
    
    write(*,*) 'Calculation is finished.'
    write(*,*) 
    write(*,*) 'Time elapsed ', delta, ' seconds.'
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
    
    do idet = 1, ndet
        fnam_temp = trim(adjustl(output_prefix)) // '_DiffPlaneTotal_Detector'
        call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
        call output_stem_image(stem_image(:,:,:,idet,:),fnam_det,probe_df)
        
        if(output_thermal) then
        fnam_temp = trim(adjustl(output_prefix)) // '_DiffPlaneElastic_Detector'
        call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
        call output_stem_image(stem_elastic_image(:,:,:,idet,:),fnam_det,probe_df)
        
        stem_inelastic_image = stem_image(:,:,:,idet,:) - stem_elastic_image(:,:,:,idet,:)
        fnam_temp = trim(adjustl(output_prefix)) // '_DiffPlaneTDS_Detector'
        call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
        call output_stem_image(stem_inelastic_image,fnam_det,probe_df)
        endif
        

    enddo

    !ionization
    if(ionization) then
        do ii=1,num_ionizations
        if(EDX) then
            filename = trim(adjustl(output_prefix)) // '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_shell_EDX'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        else
            filename = trim(adjustl(output_prefix))// '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
            stem_ion_image(:,:,:,:,ii) = stem_ion_image(:,:,:,:,ii)*eels_correction_image

            filename =  trim(adjustl(output_prefix)) //'_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS_Corrected'            
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        endif
        enddo
        if(.not.EDX) then
            filename = trim(adjustl(output_prefix)) // '_EELS_CorrectionMap' 
            call output_stem_image(eels_correction_image, filename,probe_df)
        endif
       
    endif

end subroutine qep_stem
