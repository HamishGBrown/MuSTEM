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

module cuda_setup
    
    use cudafor
    use m_user_input, only: get_input
    
    implicit none
    
    integer(8) :: usable_GPU_memory             !memory available on the GPU
    
    contains
    
    subroutine display_device_properties(i)
    
        use m_string, only: to_string
        
        implicit none
        
        integer::i
        
        integer::istat
        
        type (cudaDeviceProp) :: properties
        
         istat = cudaGetDeviceProperties(properties, i)
         
         write(*,"(' Device Number: ',i0)") i
         write(*,"('   Device name: ',a)") trim(properties%name)
         write(*,"('   Memory Clock Rate (MHz): ', i0)") properties%memoryClockRate/1000
         write(*,"('   Memory Bus Width (bits): ', i0)") properties%memoryBusWidth
         write(*,"('   Peak Memory Bandwidth (GB/s): ', f6.2)") &
           2.0*properties%memoryClockRate*(properties%memoryBusWidth/8)/10.0**6
         write(*,"('   Total Global Memory (MB): ', f10.2)") properties%totalGlobalMem/10.0**6
         
         if (properties%major.lt.2) then         
            write(*,"('   Compute capability: ', a, ' (WARNING: below minimum requirements)')") to_string(properties%major) // '.' // to_string(properties%minor)
         else
            write(*,"('   Compute capability: ', a)") to_string(properties%major) // '.' // to_string(properties%minor)
         endif
         
         


        
	end subroutine
    
    subroutine setup_GPU
    
        implicit none
    
        integer :: i, istat, nDevices
        type (cudaDeviceProp) :: properties
    
        write(*,*) '|----------------------------------|'
        write(*,*) '|           GPU selection          |'
        write(*,*) '|----------------------------------|'
        write(*,*)
    
        istat = cudaGetDeviceCount(nDevices)
    
        if (istat.ne.0) then
            write(*,*) "Your system does not appear to be CUDA-capable."
            write(*,*) "The program cannot proceed and will exit."
            pause
            stop
	    endif
    
        if (nDevices==0) then
            write(*,*) "You have no CUDA-capable devices."
            write(*,*) "The program cannot proceed and will exit."
            pause
            stop
        
	    elseif (nDevices==1) then
            write(*,*) "You have one CUDA-capable device, with the following properties:"
            write(*,*)
        
            call display_device_properties(0)
        
            istat = cudaGetDeviceProperties(properties, 0)
        
            usable_GPU_memory = properties%totalGlobalMem
        
            istat = cudaSetDevice(0)
        
            if (properties%major.lt.2) then
                write(*,*) 'Warning: this device has compute capability less than 2.0'
                write(*,*) 'which will likely result in failure.'
                write(*,*)
            endif
            
    5       write(*,*) 'Enter <0> to continue.'
            call get_input('Device used for calculation',i) 
            if (i.ne.0) goto 5
            write(*,*)
        
	    else
    
            write(*,*) "You have the following CUDA-capable devices:"
            write(*,*)
        
	        do i = 0, nDevices-1
                 call display_device_properties(i)
            enddo
        
    10      write(*,*) 'Enter the device number to use for the calculation:'
            call get_input('Device used for calculation',i) 
            write(*,*)
        
            if (i.lt.0 .or. i.ge.nDevices) goto 10
    
            istat = cudaGetDeviceProperties(properties, i)
            usable_GPU_memory = properties%totalGlobalMem
    
            istat = cudaSetDevice(i)
            
            if (properties%major.le.2) then
                write(*,*) 'Warning: this device has compute capability less than 2.0'
                write(*,*) 'which will likely result in failure.'
                write(*,*)
            endif

	    endif
    
    end subroutine
    
    
    
    subroutine GPU_memory_message(required_memory, on_the_fly)
    
        use m_precision, only: fp_kind
        
        implicit none
        
        real(fp_kind),intent(in) :: required_memory
        logical,intent(out) :: on_the_fly
        
        integer :: i_choice
        
        write(*,*) '|-----------------------|'
        write(*,*) '|      GPU memory       |'
        write(*,*) '|-----------------------|'
        write(*,*)
            
        if(required_memory/10.0**6<10.0**5) then ;write(*,'(1x, a, f9.3, a)') 'This calculation requires a GPU with roughly ', required_memory/10.0**6, ' MB of memory.'
		else; write(*,'(1x, a, f9.0, a)') 'This calculation requires a GPU with > ', 10.0**5, ' MB of memory.';endif

        write(*,*) 'The selected device has access to ', usable_GPU_memory/10.0**6, ' MB of memory.'
        write(*,*)
        
        if(required_memory.lt.usable_GPU_memory) then
            write(*,*) 'It appears that the selected GPU device satisfies this requirement,'
            write(*,*) 'and thus the calculations can proceed using precalculated potentials.'
            write(*,*) 'However, in the event that the code still crashes, the user may wish'
            write(*,*) 'to force on-the-fly calculation of potentials to see if that helps.'
            write(*,*)
            write(*,*) 'Warning: issues have been reported with on-the-fly calculations'
            write(*,*) 'and non-square arrays.'
            write(*,*)
        
10          write(*,*) '<0> Precalculated potentials'
            write(*,*) '<1> On-the-fly calculation'
            call get_input('<0> Precalculated potentials <1> On-the-fly calculation', i_choice)
            write(*,*)
            
            select case (i_choice)
                case (0)
                    on_the_fly = .false.
                    
                case (1)
                    on_the_fly = .true.
                    
                case default
                    goto 10
                    
            end select
        
        
        else
            write(*,*) 'The required memory exceeds the total available, and thus'
            write(*,*) 'calculations will proceed using on-the-fly calculated potentials.'
            write(*,*)
            write(*,*) 'Warning: issues have been reported with on-the-fly calculations'
            write(*,*) 'and non-square arrays.'
            write(*,*)
        
        
            on_the_fly = .true.
        
        endif
        
    end subroutine GPU_memory_message
    
    
    
end module
    