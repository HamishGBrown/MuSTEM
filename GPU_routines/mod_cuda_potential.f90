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

module cuda_potential
    
    use global_variables
    use cuda_array_library
    use cudafor
    
    implicit none
    
    real(fp_kind), allocatable :: ccd_slice_array(:)

    real(fp_kind), allocatable :: Volume_array(:)

    
    contains

    
    
    subroutine cuda_setup_many_phasegrate
    
        use m_precision
        use global_variables
        use m_slicing
    
        implicit none
    
        integer(4)  i,j,m,n
        integer(4)  nmax,no_atoms
        complex(fp_kind)  projected_potential(nopiy,nopix)
        complex(fp_kind)  temp(nopiy,nopix)
        real(fp_kind)  ss_temp(7)

        complex(fp_kind) :: scattering_pot(nopiy,nopix,nt)
    
        if(allocated(CCD_slice_array)) deallocate(CCD_slice_array)
        if(allocated(Volume_array)) deallocate(Volume_array)
	
        allocate(CCD_slice_array(n_slices))
        allocate(Volume_array(n_slices))
 	
        do j = 1, n_slices
            Volume_array(j) = ss_slice(7,j)
	        no_atoms = sum(nat_slice(:,j))
	     
    198	    write(6,199) no_atoms
    199     format(/,' Number of atoms in slice:',i5,/) 
		    ccd_slice_array(j) = relm / (tp * ak * ss_slice(7,j))
	    enddo
	
    end subroutine cuda_setup_many_phasegrate    
    
  attributes(global) subroutine cuda_make_atom_mask_real(tau_ss_in,input_nat_layer,atom_mask_real_d,n,m)
    
        use m_precision
    
        implicit none
    
        integer(4),value :: ix
        integer(4),value :: xpixel,ypixel,n,m
        integer(4),value :: input_nat_layer !dimensions of the tau array
        real(fp_kind) :: xpos,ypos,fracx,fracy
        real(fp_kind),device,dimension(3,input_nat_layer) :: tau_ss_in
        real(fp_kind),device,dimension(n,m) :: atom_mask_real_d
    
        integer :: istat
        
        ix = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (ix <= input_nat_layer) then
		
            !the -1 is to fix a pixel offset
            xpos = tau_ss_in(1,ix)*float(m) 
            ypos = tau_ss_in(2,ix)*float(n) 
            if(ceiling(xpos).gt.m) then
                xpos = xpos - float(m)
            elseif(floor(xpos).lt.1.0_fp_kind) then
                xpos = xpos + float(m)
            endif
            
            if(ceiling(ypos).gt.n) then
                ypos = ypos - float(n)
            elseif(floor(ypos).lt.1.0_fp_kind) then
                ypos = ypos + float(n)
            endif
            
            !fraction of the pixel top right
            xpixel = ceiling(xpos)
            ypixel = ceiling(ypos)
            fracx = mod(xpos,1.0_fp_kind)
            fracy = mod(ypos,1.0_fp_kind)
            
            if(xpixel.eq.0) xpixel = m
            if(xpixel.eq.m+1) xpixel = 1
            if(ypixel.eq.0) ypixel = n
            if(ypixel.eq.n+1) ypixel = 1
            istat = atomicAdd(atom_mask_real_d(ypixel,xpixel), fracx*fracy)
            
            !fraction of the pixel top left
            xpixel = floor(xpos)
            ypixel = ceiling(ypos)
            fracx = 1.0_fp_kind - mod(xpos,1.0_fp_kind)
            fracy = mod(ypos,1.0_fp_kind)
            
            if(xpixel.eq.0) xpixel = m
            if(xpixel.eq.m+1) xpixel = 1
            if(ypixel.eq.0) ypixel = n
            if(ypixel.eq.n+1) ypixel = 1
            istat = atomicAdd(atom_mask_real_d(ypixel,xpixel), fracx*fracy)
            
            !fraction of the pixel bottom right
            xpixel = ceiling(xpos)
            ypixel = floor(ypos)
            fracx = mod(xpos,1.0_fp_kind)
            fracy = 1.0_fp_kind - mod(ypos,1.0_fp_kind)
            
            if(xpixel.eq.0) xpixel = m
            if(xpixel.eq.m+1) xpixel = 1
            if(ypixel.eq.0) ypixel = n
            if(ypixel.eq.n+1) ypixel = 1
            istat = atomicAdd(atom_mask_real_d(ypixel,xpixel), fracx*fracy)
            
            !fraction of the pixel bottom left
            xpixel = floor(xpos)
            ypixel = floor(ypos)
            fracx = 1.0_fp_kind - mod(xpos,1.0_fp_kind)
            fracy = 1.0_fp_kind - mod(ypos,1.0_fp_kind)
            if(xpixel.eq.0) xpixel = m
            if(xpixel.eq.m+1) xpixel = 1
            if(ypixel.eq.0) ypixel = n
            if(ypixel.eq.n+1) ypixel = 1
            istat = atomicAdd(atom_mask_real_d(ypixel,xpixel), fracx*fracy)
        endif
        
end subroutine cuda_make_atom_mask_real

	attributes(device) subroutine check_pixel_wrapping(ypixel,xpixel,n,m)
			integer,intent(in)::n,m
			integer,intent(inout)::ypixel,xpixel
			xpixel = mod(xpixel-1,m)
			ypixel = mod(ypixel-1,n)
end subroutine
    
    subroutine cuda_fph_make_potential(transf_d,ccd_slice,tau_ss,nat_layer,n_sub_slice,thickness,idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
    
        use global_variables
        use m_precision
        use CUFFT
	    use cufft_wrapper
        use cudafor
        use m_slicing, only: n_slices, nat_slice, a0_slice
        use m_qep, only: displace
        use cuda_array_library, only: blocks, threads
		use output
        
        implicit none
    
        integer plan
        integer(4) :: n_sub_slice
        integer(4) :: i,j,nat_layer(nt)
        integer(4) :: idum
    
        real(fp_kind) :: tau_ss(3,nt,maxval(nat)*ifactorx*ifactory,n_slices),tau_displace(3,nt,maxval(nat_layer))
        real(fp_kind),device,allocatable :: tau_displace_d(:,:)
        real(fp_kind) :: CCD_slice,thickness
        real(fp_kind) :: image(nopiy,nopix)
        real(fp_kind) :: interaction
    
        !Device variables
        type(dim3) :: blocks_nat, threads_nat
        complex(fp_kind),device,dimension(nopiy,nopix) :: transf_d,bwl_mat_d
        complex(fp_kind),device,dimension(nopiy,nopix,nt) :: fz_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: slice_potential_d,potential_d,inverse_sinc_d,atom_mask_d
        real(fp_kind),device,dimension(nopiy,nopix) :: atom_mask_real_d, temp_d
	              
 	    do i=1,nt                                 !calculate shifted tau
            do j=1,nat_slice(i,n_sub_slice)       !loop over atoms of type i in subslice n_sub_slice
                call displace(tau_ss(1:3,i,j,n_sub_slice),tau_displace(1:3,i,j),sqrt(atf(3,i)),a0_slice,idum)
            enddo
        enddo

        slice_potential_d = 0.0_fp_kind
        
        threads_nat = dim3(1024, 1, 1)
        
        do i = 1, nt
            if(nat_layer(i).eq.0) cycle
                        
            blocks_nat = dim3(ceiling(float(nat_layer(i))/threads_nat%x), 1, 1)
            
            allocate(tau_displace_d(3,nat_layer(i)))
            tau_displace_d = tau_displace(:,i,1:nat_layer(i))
            
            atom_mask_real_d = 0.0_fp_kind
            call cuda_make_atom_mask_real<<<blocks_nat,threads_nat>>>(tau_displace_d,nat_layer(i),atom_mask_real_d,nopiy,nopix)

            call cuda_cshift_real<<<blocks,threads>>>(atom_mask_real_d, temp_d, nopiy, nopix, -1, -1)
            
            call cuda_copy_real_to_cmplx<<<blocks,threads>>>(temp_d, atom_mask_d, nopiy, nopix)           
            
            !fourier transform to convolve with scattering factor
            call cufftExec(plan,atom_mask_d,potential_d,CUFFT_FORWARD)
            
            !divide by the sinc function
            call cuda_multiplication<<<blocks,threads>>>(potential_d,inverse_sinc_d,potential_d,1.0_fp_kind,nopiy,nopix)   
        
            !multiply by the scattering factor
            call cuda_multiplication<<<blocks,threads>>>(potential_d,fz_d(:,:,i),potential_d,CCD_slice,nopiy,nopix) 
            
            !sum the complex potentials (in reciprocal space)
            call cuda_addition<<<blocks,threads>>>(slice_potential_d,potential_d,slice_potential_d,normalisation,nopiy,nopix)          
            
            deallocate(tau_displace_d)
        enddo
        
        !get realspace potential
        call cufftExec(plan,slice_potential_d,potential_d,CUFFT_INVERSE)
    
        !bandwidth limit the potential
        interaction = pi*thickness/ak1
        call cuda_complex_exponentiation<<<blocks,threads>>>(transf_d,potential_d,interaction,nopiy,nopix)  
        call cufftExec(plan,transf_d,potential_d,CUFFT_FORWARD)
        call cuda_multiplication<<<blocks,threads>>>(potential_d,bwl_mat_d,potential_d,normalisation,nopiy,nopix)  
        call cufftExec(plan,potential_d,transf_d,CUFFT_INVERSE)

    end subroutine cuda_fph_make_potential
    
    
    
    subroutine cuda_make_abs_potential(transf_d,ccd_slice,tau_ss,nat_layer,thickness,plan,fz_d,fz_dwf_d,fz_abs_d,inverse_sinc_d,bwl_mat_d,volume)
    
        use global_variables
        use m_precision
        use CUFFT
	    use cufft_wrapper
        use cudafor
        use cuda_array_library, only: blocks, threads
    
        implicit none
    
        integer plan
        integer(4) :: i,j,nat_layer(nt)
    
        real(fp_kind) :: tau_ss(3,nt,maxval(nat)*ifactorx*ifactory)
        real(fp_kind) :: CCD_slice,thickness
        real(fp_kind) :: interaction
        real(fp_kind) :: volume,V_corr
    
        !Device variables
        type(dim3) :: blocks_nat, threads_nat
        real(fp_kind),device,allocatable :: tau_d(:,:)
        complex(fp_kind),device,dimension(nopiy,nopix) :: transf_d,bwl_mat_d
        complex(fp_kind),device,dimension(nopiy,nopix,nt) :: fz_d,fz_abs_d,fz_dwf_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: slice_potential_d,potential_d,inverse_sinc_d,atom_mask_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: elastic_slice_potential_d,inelastic_slice_potential_d
        real(fp_kind),device,dimension(nopiy,nopix) :: atom_mask_real_d, temp_d
    
        complex(fp_kind),device,dimension(nopiy,nopix) :: potential_inelastic_d,potential_elastic_d
		

        elastic_slice_potential_d = 0.0_fp_kind
        inelastic_slice_potential_d = 0.0_fp_kind
        V_corr = ss(7)/volume
        
        threads_nat = dim3(1024, 1, 1)
        
        do i = 1, nt
            if(nat_layer(i).eq.0) cycle
            
            blocks_nat = dim3(ceiling(float(nat_layer(i))/threads_nat%x), 1, 1)
            
            allocate(tau_d(3,nat_layer(i)))
            tau_d = tau_ss(:,i,1:nat_layer(i))

            atom_mask_real_d = 0.0_fp_kind
            call cuda_make_atom_mask_real<<<blocks_nat,threads_nat>>>(tau_d,nat_layer(i),atom_mask_real_d,nopiy,nopix)
            call cuda_cshift_real<<<blocks,threads>>>(atom_mask_real_d, temp_d, nopiy, nopix, -1, -1)
            call cuda_copy_real_to_cmplx<<<blocks,threads>>>(temp_d, atom_mask_d, nopiy, nopix)           
            
            deallocate(tau_d)
            
            !fourier transform to convolve with scattering factor
            call cufftExec(plan,atom_mask_d,potential_d,CUFFT_FORWARD)
            
            !multiply by the scattering factor
            !Elastic potential
            call cuda_multiplication<<<blocks,threads>>>(potential_d,fz_d(:,:,i),potential_elastic_d,1.0_fp_kind,nopiy,nopix)  
            
            !DWF potential
            call cuda_multiplication<<<blocks,threads>>>(potential_elastic_d,fz_dwf_d(:,:,i),potential_elastic_d,CCD_slice,nopiy,nopix) 
            
            !absorptive potential
            call cuda_multiplication<<<blocks,threads>>>(potential_d,fz_abs_d(:,:,i),potential_inelastic_d,V_corr,nopiy,nopix)
            
            !sum the complex potentials (in reciprocal space)
            call cuda_addition<<<blocks,threads>>>(elastic_slice_potential_d,potential_elastic_d,elastic_slice_potential_d,normalisation,nopiy,nopix) 
            call cuda_addition<<<blocks,threads>>>(inelastic_slice_potential_d,potential_inelastic_d,inelastic_slice_potential_d,normalisation,nopiy,nopix) 
        enddo
        
        !divide by the sinc function
        call cuda_addition<<<blocks,threads>>>(inelastic_slice_potential_d,elastic_slice_potential_d,slice_potential_d,1.0_fp_kind,nopiy,nopix)
        call cuda_multiplication<<<blocks,threads>>>(slice_potential_d,inverse_sinc_d,slice_potential_d,1.0_fp_kind,nopiy,nopix) 
        
        !get realspace potential
        call cufftExec(plan,slice_potential_d,potential_d,CUFFT_INVERSE)
    
        !bandwidth limit the potential
        interaction = pi*thickness/ak1
        call cuda_complex_exponentiation<<<blocks,threads>>>(transf_d,potential_d,interaction,nopiy,nopix)  
        call cufftExec(plan,transf_d,potential_d,CUFFT_FORWARD)
        call cuda_multiplication<<<blocks,threads>>>(potential_d,bwl_mat_d,potential_d,normalisation,nopiy,nopix)  
        call cufftExec(plan,potential_d,transf_d,CUFFT_INVERSE)

    end subroutine cuda_make_abs_potential
    
    

    subroutine cuda_make_adf_potential(real_inelastic_slice_potential_d,tau_ss,nat_layer,plan,fz_adf_d,inverse_sinc_d,volume)
    
        use global_variables
        use m_precision
        use CUFFT
	    use cufft_wrapper
        use cudafor
        use cuda_array_library, only: blocks, threads
        
        implicit none
    
        integer plan
        integer(4) :: i,j,nat_layer(nt)
    
        real(fp_kind) :: tau_ss(3,nt,maxval(nat)*ifactorx*ifactory)
        real(fp_kind) :: volume,V_corr
    
        !Device variables
        type(dim3) :: blocks_nat, threads_nat
        real(fp_kind),device,allocatable :: tau_d(:,:)
        real(fp_kind),device,dimension(nopiy,nopix) :: real_inelastic_slice_potential_d
        complex(fp_kind),device,dimension(nopiy,nopix,nt) :: fz_adf_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: slice_potential_d,potential_d,inverse_sinc_d,atom_mask_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: inelastic_slice_potential_d
        complex(fp_kind),dimension(nopiy,nopix) :: atom_mask
        real(fp_kind),device,dimension(nopiy,nopix) :: atom_mask_real_d, temp_d
    
        V_corr = ss(7)/volume
        inelastic_slice_potential_d = 0.0_fp_kind
        
        threads_nat = dim3(1024, 1, 1)
        
        do i =1,nt
            if(nat_layer(i).eq.0) cycle
            
            blocks_nat = dim3(ceiling(float(nat_layer(i))/threads_nat%x), 1, 1)
            
            allocate(tau_d(3,nat_layer(i)))
            tau_d = tau_ss(:,i,1:nat_layer(i))
            
            atom_mask_real_d = 0.0_fp_kind
            call cuda_make_atom_mask_real<<<blocks_nat,threads_nat>>>(tau_d,nat_layer(i),atom_mask_real_d,nopiy,nopix)
            call cuda_cshift_real<<<blocks,threads>>>(atom_mask_real_d, temp_d, nopiy, nopix, -1, -1)
            call cuda_copy_real_to_cmplx<<<blocks,threads>>>(temp_d, atom_mask_d, nopiy, nopix)           
            
            deallocate(tau_d)
            !fourier transform to convolve with scattering factor
            call cufftExec(plan,atom_mask_d,potential_d,CUFFT_FORWARD)
            !multiply by the scattering factor
            !Elastic potential
            call cuda_multiplication<<<blocks,threads>>>(potential_d,fz_adf_d(:,:,i),potential_d,V_corr,nopiy,nopix)  
            call cuda_addition<<<blocks,threads>>>(inelastic_slice_potential_d,potential_d,inelastic_slice_potential_d,normalisation,nopiy,nopix) 
        enddo
        !divide by the sinc function
        call cuda_multiplication<<<blocks,threads>>>(inelastic_slice_potential_d,inverse_sinc_d,potential_d,1.0_fp_kind,nopiy,nopix)   
        !get realspace potential
        call cufftExec(plan,potential_d,inelastic_slice_potential_d,CUFFT_INVERSE)
        call cuda_take_real<<<blocks,threads>>>(inelastic_slice_potential_d,real_inelastic_slice_potential_d,nopiy,nopix)   
    
    end subroutine cuda_make_adf_potential
    
    
    
    subroutine cuda_make_ion_potential(real_inelastic_slice_potential_d,tau_ss,nat_layer,plan,fz_mu_d,inverse_sinc_d,volume)
    
        use global_variables
        use m_precision
        use CUFFT
	    use cufft_wrapper
        use cudafor
        use cuda_array_library, only: blocks, threads
        
        implicit none
    
        integer plan
        integer(4) :: j,nat_layer,kval
    
        real(fp_kind) :: tau_ss(3,maxval(nat)*ifactorx*ifactory)
        real(fp_kind) :: volume,absorptive_scale
    
        !Device variables
        type(dim3) :: blocks_nat, threads_nat
		!max(nat_layer,1) avoids a bug where the program attempts to allocate zero memory
        real(fp_kind),device :: tau_d(3,max(nat_layer,1))
        real(fp_kind),device,dimension(nopiy,nopix) :: real_inelastic_slice_potential_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: fz_mu_d,fz_dwf_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: slice_potential_d,potential_d,inverse_sinc_d,atom_mask_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: inelastic_slice_potential_d
        complex(fp_kind),dimension(nopiy,nopix) :: atom_mask
        real(fp_kind),device,dimension(nopiy,nopix) :: atom_mask_real_d, temp_d
    
        absorptive_scale = 1.0_fp_kind/volume
        inelastic_slice_potential_d = 0.0_fp_kind
		if(nat_layer>1) then
        
        !make the atom mask
        threads_nat = dim3(1024, 1, 1)
        blocks_nat = dim3(ceiling(float(nat_layer)/threads_nat%x), 1, 1)
        
        tau_d = tau_ss(:,1:nat_layer)
        
        atom_mask_real_d = 0.0_fp_kind
        call cuda_make_atom_mask_real<<<blocks_nat,threads_nat>>>(tau_d,nat_layer,atom_mask_real_d,nopiy,nopix)
        call cuda_cshift_real<<<blocks,threads>>>(atom_mask_real_d, temp_d, nopiy, nopix, -1, -1)
        call cuda_copy_real_to_cmplx<<<blocks,threads>>>(temp_d, atom_mask_d, nopiy, nopix)           
        
        !fourier transform to convolve with scattering factor
        call cufftExec(plan,atom_mask_d,potential_d,CUFFT_FORWARD)
                
        !Elastic potential
        call cuda_multiplication<<<blocks,threads>>>(potential_d,fz_mu_d,potential_d,1.0_fp_kind,nopiy,nopix)  
        
        !DWF potential
        !call cuda_multiplication<<<blocks,threads>>>(potential_d,fz_dwf_d,potential_d,absorptive_scale,nopiy,nopix) 
        call cuda_addition<<<blocks,threads>>>(inelastic_slice_potential_d,potential_d,inelastic_slice_potential_d,normalisation*absorptive_scale,nopiy,nopix) 
        
        !divide by the sinc function
        call cuda_multiplication<<<blocks,threads>>>(inelastic_slice_potential_d,inverse_sinc_d,potential_d,1.0_fp_kind,nopiy,nopix)   
        
        !get realspace potential
        call cufftExec(plan,potential_d,inelastic_slice_potential_d,CUFFT_INVERSE)

		endif
        call cuda_take_real<<<blocks,threads>>>(inelastic_slice_potential_d,real_inelastic_slice_potential_d,nopiy,nopix)   
    
    end subroutine cuda_make_ion_potential
    
    

    attributes(global) subroutine cuda_site_factor(site_factor, tau, g_vec_array, n, m)
    
        implicit none
     
        integer(4), value :: n, m, ix, iy, i
        complex(fp_kind),dimension(n,m) :: site_factor
        integer,dimension(3,n,m) :: g_vec_array
        real(fp_kind) :: tau(:,:)
        real(fp_kind),parameter :: tp = atan(1.0_fp_kind)*8.0_fp_kind 
    
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (ix <= m) then
             if (iy <= n) then
                site_factor(iy,ix) = 0.0_fp_kind
            
                do i = 1, size(tau, 2)
                    site_factor(iy,ix) = site_factor(iy,ix) + exp(cmplx(0.0_fp_kind, -tp*(g_vec_array(1,iy,ix)*tau(1,i)+g_vec_array(2,iy,ix)*tau(2,i)), fp_kind))
                enddo
            endif
        endif   
    
    end subroutine cuda_site_factor

    
    
    end module
    