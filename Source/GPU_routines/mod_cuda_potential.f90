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
		use m_multislice,only:n_slices,nat_slice,ss_slice
    
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
! 
!     attributes(device) subroutine atAddComplex( a, b)
!
!	 implicit none
!
!	 complex(fp_kind) :: a
!	 complex(fp_kind) :: b
!
!     
!      call atomicAdd(real(a), real(b))
!      call atomicAdd(imag(a), imag(b))
!    end subroutine 
!
!   attributes(global) subroutine cuda_make_atom_mask_complex(tau_ss_in,input_nat_layer,atom_mask_d,n,m)
!    
!        use m_precision
!    
!        implicit none
!    
!        integer(4),value :: ix
!        integer(4),value :: xpixel,ypixel,n,m
!        integer(4),value :: input_nat_layer !dimensions of the tau array
!        real(fp_kind) :: xpos,ypos,fracx,fracy
!        real(fp_kind),device,dimension(3,input_nat_layer) :: tau_ss_in
!        complex(fp_kind),device,dimension(n,m) :: atom_mask_d
!    
!        integer :: istat
!        
!        ix = (blockIdx%x-1)*blockDim%x + threadIdx%x
!    
!        if (ix <= input_nat_layer) then
!		
!            !the -1 is to fix a pixel offset
!            xpos = tau_ss_in(1,ix)*float(m)
!            ypos = tau_ss_in(2,ix)*float(n)
!            if(ceiling(xpos).gt.m) then
!                xpos = xpos - float(m)
!            elseif(floor(xpos).lt.1.0_fp_kind) then
!                xpos = xpos + float(m)
!            endif
!            
!            if(ceiling(ypos).gt.n) then
!                ypos = ypos - float(n)
!            elseif(floor(ypos).lt.1.0_fp_kind) then
!                ypos = ypos + float(n)
!            endif
!            
!            !fraction of the pixel top right
!            xpixel = ceiling(xpos)
!            ypixel = ceiling(ypos)
!            fracx = mod(xpos,1.0_fp_kind)
!            fracy = mod(ypos,1.0_fp_kind)
!            
!            if(xpixel.eq.0) xpixel = m
!            if(xpixel.eq.m+1) xpixel = 1
!            if(ypixel.eq.0) ypixel = n
!            if(ypixel.eq.n+1) ypixel = 1
!            call atAddComplex(atom_mask_d(ypixel,xpixel), fracx*fracy)
!            
!            !fraction of the pixel top left
!            xpixel = floor(xpos)
!            ypixel = ceiling(ypos)
!            fracx = 1.0_fp_kind - mod(xpos,1.0_fp_kind)
!            fracy = mod(ypos,1.0_fp_kind)
!            
!            if(xpixel.eq.0) xpixel = m
!            if(xpixel.eq.m+1) xpixel = 1
!            if(ypixel.eq.0) ypixel = n
!            if(ypixel.eq.n+1) ypixel = 1
!            call atAddComplex(atom_mask_d(ypixel,xpixel), fracx*fracy)
!            
!            !fraction of the pixel bottom right
!            xpixel = ceiling(xpos)
!            ypixel = floor(ypos)
!            fracx = mod(xpos,1.0_fp_kind)
!            fracy = 1.0_fp_kind - mod(ypos,1.0_fp_kind)
!            
!            if(xpixel.eq.0) xpixel = m
!            if(xpixel.eq.m+1) xpixel = 1
!            if(ypixel.eq.0) ypixel = n
!            if(ypixel.eq.n+1) ypixel = 1
!            call atAddComplex(atom_mask_d(ypixel,xpixel), fracx*fracy)
!            
!            !fraction of the pixel bottom left
!            xpixel = floor(xpos)
!            ypixel = floor(ypos)
!            fracx = 1.0_fp_kind - mod(xpos,1.0_fp_kind)
!            fracy = 1.0_fp_kind - mod(ypos,1.0_fp_kind)
!            if(xpixel.eq.0) xpixel = m
!            if(xpixel.eq.m+1) xpixel = 1
!            if(ypixel.eq.0) ypixel = n
!            if(ypixel.eq.n+1) ypixel = 1
!            call atAddComplex(atom_mask_d(ypixel,xpixel), fracx*fracy)
!        endif
!        
!end subroutine cuda_make_atom_mask_complex
 
    
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
            xpos = tau_ss_in(1,ix)*float(m)+1.0_fp_kind
            ypos = tau_ss_in(2,ix)*float(n)+1.0_fp_kind
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
        use m_numerical_tools, only: displace
        use cuda_array_library, only: blocks, threads
		use output
		use m_multislice,only:n_slices,nat_slice,a0_slice
        
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
            call cuda_copy_real_to_cmplx<<<blocks,threads>>>(atom_mask_real_d, potential_d, nopiy, nopix)           
            
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
    
  
    subroutine cuda_on_the_fly_abs_multislice(psi_d,ccd_slice,tau_ss,nat_layer,thickness,even_slicing,plan,Vg,volume,propy_d, propx_d,prop_d)
    
        use global_variables
        use m_precision
        use CUFFT
	    use cufft_wrapper
        use cudafor
        use cuda_array_library, only: blocks, threads
		use output
    
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
		complex(fp_kind),device,dimension(nopiy,nopix),intent(inout) :: psi_d
        complex(fp_kind),device,dimension(nopiy,nopix) :: transf_d,potential_d,Vg_d
		complex(fp_kind),device,intent(in),optional::propy_d(nopiy),propx_d(nopix),prop_d(nopiy,nopix)
		complex(fp_kind),device,allocatable::propp_d(:,:)
        complex(fp_kind),dimension(nopiy,nopix,nt),intent(in) :: Vg
        real(fp_kind),device,dimension(nopiy,nopix) :: atom_mask_real_d
		logical,intent(in)::even_slicing
		complex(fp_kind),dimension(nopiy,nopix)::out

		transf_d = 0
        V_corr = ss(7)/volume
        
        threads_nat = dim3(1024, 1, 1)
        
        do i = 1, nt
            if(nat_layer(i).eq.0) cycle
            
            blocks_nat = dim3(ceiling(float(nat_layer(i))/threads_nat%x), 1, 1)
            
            allocate(tau_d(3,nat_layer(i)))
            tau_d = tau_ss(:,i,1:nat_layer(i))

            atom_mask_real_d = 0.0_fp_kind
            call cuda_make_atom_mask_real<<<blocks_nat,threads_nat>>>(tau_d,nat_layer(i),atom_mask_real_d,nopiy,nopix)

            call cuda_copy_real_to_cmplx<<<blocks,threads>>>(atom_mask_real_d, potential_d, nopiy, nopix)           
            deallocate(tau_d)

            !fourier transform to convolve with scattering factor
            call cufftExec(plan,potential_d,potential_d,CUFFT_FORWARD)

            !multiply by the scattering factor
            !Elastic potential
			Vg_d = Vg(:,:,i)
            call cuda_multiplication<<<blocks,threads>>>(potential_d,Vg_d,potential_d,CCD_slice,nopiy,nopix)
            
            !sum the complex potentials (in reciprocal space)
            call cuda_addition<<<blocks,threads>>>(transf_d,potential_d,transf_d,normalisation,nopiy,nopix)
			
			            
        enddo
        
        !get realspace potential
        call cufftExec(plan,transf_d,transf_d,CUFFT_INVERSE)

        !bandwidth limit the potential
        interaction = pi*thickness/ak1
        call cuda_complex_exponentiation<<<blocks,threads>>>(transf_d,transf_d,interaction,nopiy,nopix)  
        call cufftExec(plan,transf_d,transf_d,CUFFT_FORWARD)
		call cuda_apply_Band_width_limit<<<blocks,threads>>>(transf_d,nopiy,nopix,normalisation)
		
		call cufftExec(plan,transf_d,transf_d,CUFFT_INVERSE)
		!Multiply transmission function by psi
        
		call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d, psi_d,1.0_fp_kind,nopiy,nopix)
		
		!Apply free space propagator
		call cufftExec(plan,psi_d,psi_d,CUFFT_FORWARD)
        if(present(propy_d).and.present(propx_d)) then
			call cuda_multiplication<<<blocks,threads>>>(psi_d,propy_d,propx_d, psi_d,normalisation,nopiy,nopix)
		elseif(present(prop_d)) then
			if(even_slicing) then
				call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d, psi_d,normalisation,nopiy,nopix)
			else
				allocate(propp_d(nopiy,nopix))
				call cuda_complex_exponentiation<<<blocks,threads>>>(propp_d, prop_d,thickness, nopiy, nopix)
				call cuda_apply_Band_width_limit<<<blocks,threads>>>(propp_d, nopiy,nopix,1.0_fp_kind)
				call cuda_multiplication<<<blocks,threads>>>(psi_d,propp_d, psi_d,normalisation,nopiy,nopix)
				
			endif
		else
			stop
		endif
			
        call cufftExec(plan,psi_d,psi_d,CUFFT_INVERSE)

		

    end subroutine cuda_on_the_fly_abs_multislice   

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
    
	function cuda_adf_crossection_on_the_fly(psi_d,tau_ss,nat_layer,plan,fz_adf,volume,prop_distance)
    
        use global_variables
        use m_precision
        use CUFFT
	    use cufft_wrapper
        use cudafor
        use cuda_array_library, only: blocks, threads
		use output
        
        implicit none
    
        integer plan
        integer(4) :: i,j,nat_layer(nt)
    
        real(fp_kind) :: tau_ss(3,nt,maxval(nat)*ifactorx*ifactory)
        real(fp_kind) :: volume,prop_distance,V_corr
    
        !Device variables
        type(dim3) :: blocks_nat, threads_nat
        real(fp_kind),device,allocatable :: tau_d(:,:)
        real(fp_kind),device,dimension(nopiy,nopix) :: real_inelastic_slice_potential_d
		complex(fp_kind),dimension(nopiy,nopix,nt),intent(in):: fz_adf
        complex(fp_kind),device,dimension(nopiy,nopix) :: potential_d,fz_adf_d
		complex(fp_kind),dimension(nopiy,nopix),intent(in),device::psi_d
		complex(fp_kind)::out(nopiy,nopix)
		
		real(fp_kind)::cuda_adf_crossection_on_the_fly
    
        V_corr = ss(7)/volume
        cuda_adf_crossection_on_the_fly = 0.0_fp_kind
        threads_nat = dim3(1024, 1, 1)
        
        do i =1,nt
            if(nat_layer(i).eq.0) cycle
            
            blocks_nat = dim3(ceiling(float(nat_layer(i))/threads_nat%x), 1, 1)
            
            allocate(tau_d(3,nat_layer(i)))
            tau_d = tau_ss(:,i,1:nat_layer(i))
            
            real_inelastic_slice_potential_d = 0.0_fp_kind
            call cuda_make_atom_mask_real<<<blocks_nat,threads_nat>>>(tau_d,nat_layer(i),real_inelastic_slice_potential_d,nopiy,nopix)
            call cuda_copy_real_to_cmplx<<<blocks,threads>>>(real_inelastic_slice_potential_d, potential_d, nopiy, nopix)           
            deallocate(tau_d)
			
            !fourier transform to convolve with scattering factor
            call cufftExec(plan,potential_d,potential_d,CUFFT_FORWARD)
            !multiply by the scattering factor
            !Elastic potential
			fz_adf_d = fz_adf(:,:,i)
            call cuda_multiplication<<<blocks,threads>>>(potential_d,fz_adf_d,potential_d,V_corr*normalisation*prop_distance,nopiy,nopix)  
			call cufftExec(plan,potential_d,potential_d,CUFFT_INVERSE)
			
			call cuda_take_real<<<blocks,threads>>>(potential_d,real_inelastic_slice_potential_d,nopiy,nopix)
			
			cuda_adf_crossection_on_the_fly = cuda_adf_crossection_on_the_fly + local_potential_overlap(psi_d,real_inelastic_slice_potential_d,nopiy,nopix)
        enddo
        !get realspace potential
    
    end function cuda_adf_crossection_on_the_fly
    
	attributes(host) function local_potential_overlap(psi_d,pot,nopiy,nopix)
		
		implicit none
    
        integer(4) l,m,nopiy,nopix
        real(fp_kind) :: local_potential_overlap
        real(fp_kind),device :: local_potential_overlap_d
        complex(fp_kind), device, dimension(nopiy,nopix),intent(in) :: psi_d
		real(fp_kind),device,dimension(nopiy,nopix)::pot
    
        local_potential_overlap_d = 0.0_fp_kind
    
        !$cuf kernel do (2) <<<*,*>>>
        do m = 1, nopix
            do l = 1, nopiy
                local_potential_overlap_d = local_potential_overlap_d + abs(psi_d(l,m))**2*pot(l,m)
            enddo
        enddo
        
        local_potential_overlap = local_potential_overlap_d
	
	end function
    
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

	subroutine cuda_band_width_limit(arrayin,nopiy,nopix,szey,szex,normalisation,fftin,fftout,plan)
	
        use CUFFT
	    use cufft_wrapper
	complex(fp_kind),intent(inout),device::arrayin(nopiy,nopix)
	integer*4,intent(in)::nopiy,nopix
	integer,intent(in)::plan
    real(fp_kind),intent(in) :: szey,szex,normalisation
	logical,intent(in)::fftin,fftout

	call cufftExec(plan,arrayin,arrayin,CUFFT_FORWARD)
        !call cuda_multiplication<<<blocks,threads>>>(arrayin,bwl_mat_d,transf_d,normalisation,nopiy,nopix)  
        call cufftExec(plan,arrayin,arrayin,CUFFT_INVERSE)
    end subroutine
    
	attributes(global) subroutine cuda_apply_Band_width_limit(arrayin, n,m,scale)
    
        implicit none
     
        complex(fp_kind),intent(inout),device::arrayin(:,:)
		integer*4,intent(in), value::n,m
		real(fp_kind),value,intent(in):: scale
		
		integer(4), value :: ix,iy
		real(fp_kind),value :: iix,iiy
    
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (ix <= m) then
             if (iy <= n) then
                iix = abs(modulo(ix+m/2,m)-m/2)/(m/2.0_fp_kind)
				iiy = abs(modulo(iy+n/2,n)-n/2)/(n/2.0_fp_kind)
				iix = sqrt(iix**2+iiy**2)
				if(.not.(iix<2.0_fp_kind/3)) then
					arrayin(iy,ix) = 0
				else
					arrayin(iy,ix) = arrayin(iy,ix)*scale
				endif
            endif
        endif   
    
    end subroutine cuda_apply_Band_width_limit

    end module
    