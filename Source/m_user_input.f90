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

module m_user_input

    use m_precision

    implicit none

    interface get_input
        module procedure get_input_int
        module procedure get_input_two_ints
        module procedure get_input_int_array
        module procedure get_input_char
        module procedure get_input_real
        module procedure get_input_two_reals
        module procedure get_input_real_array
    end interface
    
    interface write_to_in_file
        module procedure write_to_in_file_int
        module procedure write_to_in_file_two_ints
        module procedure write_to_in_file_int_array
        module procedure write_to_in_file_char
        module procedure write_to_in_file_real
        module procedure write_to_in_file_two_reals
        module procedure write_to_in_file_real_array
    end interface



    integer(4) :: output_file_number

    
    
    integer(4) :: input_file_number, line_no

    
    integer(4),parameter :: control_file_number = 90
    integer(4),parameter :: in_file_number = 91
    
    contains

    subroutine make_new_user_input_file()
    implicit none
    integer*4::ichoice,ichoice2
    character(120)::temp_string
    logical::exists,exists2
    
    open(unit=control_file_number, file="user_input.txt", status='new')

            ichoice=-1
            do while(ichoice<1.or.ichoice>3)
                
                write(*,*) 'Please select from the following options:'
                write(*,*) '<1> Use the program in interactive mode'
                write(*,*) '<2> Record a new driving file'
                write(*,*) '<3> Play an already existing driving file'
                read(*,*) ichoice
            
            
            if(ichoice==1) write(control_file_number,*) 'interactive'
            
            if(ichoice==2) then
                ichoice2 = 2
                do while(ichoice2==2)
                    write(*,*) 'Please input the name of the new driving file'
                    read(*,*) temp_string
                    inquire(file=temp_string,exist = exists2)
                    if (exists2) then
                        write(*,*) trim(adjustl(temp_string)),' already exists.'
                        ichoice2=-1
                        do while(ichoice2<1.or.ichoice2>2)
                            write(*,*) 'Input <1> to overwrite existing driving file' 
                            write(*,*) 'Input <2> to choose a new driving file name'
                            read(*,*) ichoice2
                        enddo
                        if(ichoice2==1) write(control_file_number,*) 'record overwrite',char(10),trim(adjustl(temp_string))
                    elseif(.not.exists2) then
                        write(control_file_number,*) 'record',char(10),trim(adjustl(temp_string))
                        ichoice2=1
                    endif
                enddo
            elseif(ichoice==3) then
                ichoice2=1
                do while(ichoice2==1)
                    write(*,*) 'Please input the name of the driving file to play'
                    read(*,*) temp_string
                    inquire(file=temp_string,exist = exists2)
                    if(exists2) then
                        write(control_file_number,*) 'play',char(10),trim(adjustl(temp_string))
                        ichoice2=2
                    else
                        do while(ichoice2<1.or.ichoice2>2)
                            write(*,*) "Couldn't find ",trim(adjustl(temp_string))," options:"
                            write(*,*) "<1> Enter a new driver file name"
                            write(*,*) "<2> Return to the previous menu (eg. to record a file of that name)"
                            read(*,*) ichoice2
                        enddo
                        if(ichoice2==2) ichoice=-1
                    endif
                enddo
            endif
            enddo
            close(unit=control_file_number)
            end subroutine
            
    function init_input()
        use m_string,only:to_lower
        implicit none
        
        integer*4::init_input,io,i,ichoice,ichoice2,reason
        character(120) :: temp_string
        logical::exists,exists2
        line_no = 0
        init_input = 1
        
        inquire(file="user_input.txt",exist = exists)
        !If the user_input file does not already exist then make new one
        if(.not.exists) then 
            write(*,*);write(*,*) 'No "user_input.txt" file found, one will be created.'
            call make_new_user_input_file()
        endif
        
        open(unit=control_file_number, file="user_input.txt", status='old', err = 997)

        read(control_file_number, '(A120)',IOSTAT=Reason) temp_string
        
        if(reason<0) then
            write(*,*) char(10),'"user_input.txt" seems to be empty, creating a new one'
            close(control_file_number, status='delete')
            call make_new_user_input_file()
            open(unit=control_file_number, file="user_input.txt", status='old', err = 997)
            read(control_file_number, '(A120)',IOSTAT=Reason) temp_string
        endif
            

        if(to_lower(trim(adjustl(temp_string))) .eq. "interactive") then
            input_file_number = 5
            call init_in_file(-1)

        elseif(to_lower(trim(adjustl(temp_string))) .eq. "record") then
            input_file_number = 5
            read(control_file_number, '(A120)') temp_string
            inquire(file=trim(temp_string),exist = exists)
            if(exists) then
                write(*,*) 'File "',trim(adjustl(temp_string)),'", listed in "user_input.txt"'," already exists, and you've"
                write(*,*) 'attempted to record a file of the same name, if you would like to'
                write(*,*) 'overwrite the already existing file change "record" to "record overwrite"'
                write(*,*) 'in user_input.txt and rerun muSTEM.'
                pause
                stop
            endif
            open(unit=in_file_number, file=trim(temp_string), status='new')
            call init_in_file(in_file_number)

        elseif(to_lower(trim(adjustl(temp_string))) .eq. "record overwrite") then
          input_file_number = 5
          read(control_file_number, '(A120)') temp_string
          open(unit=in_file_number, file=trim(temp_string), status='replace')
          call init_in_file(in_file_number)

        elseif(to_lower(trim(adjustl(temp_string))) .eq. "play") then
          input_file_number = in_file_number
          read(control_file_number, '(A120)') temp_string
          !open(unit=in_file_number, file=trim(temp_string), status='old', err = 998)
          call init_in_file(-1)
          
          
        elseif(to_lower(trim(adjustl(temp_string))).eq. "play all") then
          input_file_number = in_file_number
          init_input = 0
          do
            read(control_file_number, '(A120)',iostat = io) temp_string
            if(io==0) init_input = init_input + 1
            if(io<0) exit
          enddo
          call init_in_file(-1)
        else
            write(*,*) "ERROR: The first line of user_input.txt must be one of"
            write(*,*) "    interactive"
            write(*,*) "    record"
            write(*,*) "    record overwrite"
            write(*,*) "    play"
            write(*,*) "    play all"
            write(*,*) 
            
            pause
            stop
            
        endif
        
        close(control_file_number)
        return
        
997     write(*,*) 'ERROR: USER_INPUT.TXT is missing '
        goto 999
998     write(*,*) 'ERROR:',trim(temp_string),' does not exist '
999     write(*,*) 'Defaulting to RECORD OVERWRITE status'  
        write(*,*) 'Enter the filename to record the simulation options:'
        read(*,*) temp_string
        input_file_number = 5
        open(unit=in_file_number, file=trim(temp_string), status='replace')
        call init_in_file(in_file_number)

    end function
    
    function get_driver_file(ifile)
        integer*4,intent(in)::ifile
        integer*4::i
        character(120)::get_driver_file
        line_no = 0
        open(unit=control_file_number, file="user_input.txt", status='old', err = 997)
        do i=1,ifile+1
            read(control_file_number, '(A120)') get_driver_file
        enddo
        close(control_file_number)
        return
997     stop 'ERROR: USER_INPUT.TXT is missing '
    
    end function
    
	function get_string_from_file(input_filenumber_,line_no,formatter)

		character(128)::get_string_from_file
		character*(*),optional,intent(in)::formatter
		integer*4,intent(in)::input_filenumber_
		integer*4,intent(inout)::line_no

		logical::commented
		integer*4::excl,iostat

		commented = .true.

	    do while(commented)
			if(present(formatter)) then
				read(input_filenumber_, formatter) get_string_from_file	
			else
				read(input_filenumber_, '(a128)') get_string_from_file
			endif
			!Check for comment character and cycle through file if lines are commented
			excl = index(adjustl(get_string_from_file),'!')
			commented = excl==1
			line_no = line_no + 1
		enddo
		
		if (excl>0) get_string_from_file = get_string_from_file(:excl-1)

	end function
	
    subroutine get_input_int(prompt, num)

        implicit none

        character(*) :: prompt
        integer(4) :: num

        character(128)::s
        integer :: iostat
        
        call test_prompt(prompt)
               
5       write(*,'(1x, a)', advance='no') '> '
        s = get_string_from_file(input_file_number,line_no)
        
        if (input_file_number.eq.in_file_number) write(*,*) trim(adjustl(s))

        if (len_trim(s).eq.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        s = adjustl(remove_tabs(s))

        
        read(s, '(i10)', iostat=iostat) num
        
        if (iostat.ne.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        !line_no = line_no + 1
        call write_to_in_file(prompt, num)

    end subroutine



    subroutine get_input_char(prompt, s)

        implicit none

        character(*) :: prompt
        character(*) :: s  

        call test_prompt(prompt)
          
5       write(*,'(1x, a)', advance='no') '> '
        s = get_string_from_file(input_file_number,line_no)
        
        if (input_file_number.eq.in_file_number) write(*,*) trim(adjustl(s))

        !line_no = line_no + 1
        call write_to_in_file(prompt, s)

    end subroutine



    subroutine get_input_real(prompt, num)

        implicit none

        character(*) :: prompt
        real(fp_kind) :: num

        character(128) :: s
        integer :: iostat

        call test_prompt(prompt)
        
5       write(*,'(1x, a)', advance='no') '> '
        s = get_string_from_file(input_file_number,line_no)
        
        if (input_file_number.eq.in_file_number) write(*,*) trim(adjustl(s))

        if (len_trim(s).eq.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        if (trim(adjustl(s)).eq.'f' .or. trim(adjustl(s)).eq.'t') then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        s = adjustl(remove_tabs(s))
        
        read(s, *, iostat=iostat) num
        
        if (iostat.ne.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        !line_no = line_no + 1
        call write_to_in_file(prompt, num)

    end subroutine



    subroutine get_input_two_reals(prompt, num1, num2)

        implicit none

        character(*) :: prompt
        real(fp_kind) num1, num2

        character(128)::s
        integer :: iostat
        
        call test_prompt(prompt)
        
5       write(*,'(1x, a)', advance='no') '> '
        s = get_string_from_file(input_file_number,line_no)
        
        if (input_file_number.eq.in_file_number) write(*,*) trim(adjustl(s))

        if (len_trim(s).eq.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        if (index(s, 'f').gt.0 .or. index(s, 'g').gt.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        s = adjustl(remove_tabs(s))
        
        read(s, *, iostat=iostat) num1, num2
        
        if (iostat.ne.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        !line_no = line_no + 1
        call write_to_in_file(prompt, num1, num2)

    end subroutine


      
    subroutine get_input_int_array(prompt, array, length)

        implicit none
        
        character(*) :: prompt
        integer(4) :: length
        integer(4) :: array(length)
        
        character(1024) :: s
        integer :: iostat
        integer(4) :: array2(length+1)
        
        call test_prompt(prompt)
        
5       write(*,'(1x, a)', advance='no') '> '
		s = get_string_from_file(input_file_number,line_no,formatter = '(a128)')
        
        if (input_file_number.eq.in_file_number) write(*,*) trim(adjustl(s))

        if (len_trim(s).eq.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        if (index(s, 'f').gt.0 .or. index(s, 'g').gt.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        s = adjustl(remove_tabs(s))
        
        !read(s, *, iostat=iostat) array2
        !
        !if (iostat.eq.0) then
        !    write(*,*) 'Invalid input, try again.'
        !    goto 5
        !endif
        
        read(s, *, iostat=iostat) array
        
        if (iostat.ne.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        !line_no = line_no + 1
        call write_to_in_file(prompt, array, length)

    end subroutine



    subroutine get_input_two_ints(prompt, num1, num2)

        implicit none

        character(*) :: prompt
        integer(4) num1, num2

        character(128)::s
        integer :: iostat
        
        call test_prompt(prompt)
        
5       write(*,'(1x, a)', advance='no') '> '
        s = get_string_from_file(input_file_number,line_no)
        
        if (input_file_number.eq.in_file_number) write(*,*) trim(adjustl(s))

        if (len_trim(s).eq.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        if (index(s, 'f').gt.0 .or. index(s, 'g').gt.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        s = adjustl(remove_tabs(s))
        
        read(s, *, iostat=iostat) num1, num2
        
        if (iostat.ne.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        !line_no = line_no + 1
        call write_to_in_file(prompt, num1, num2)

    end subroutine
      
      
      
    subroutine get_input_real_array(prompt, array, length)
    
        implicit none
        
        character(*) :: prompt
        integer(4) :: length
        real(fp_kind) :: array(length)
        
        character(1024) :: s
        integer :: iostat
        real(fp_kind) :: array2(length+1)
        
        call test_prompt(prompt)
        
5       write(*,'(1x, a)', advance='no') '> '
        s = get_string_from_file(input_file_number,line_no)
        
        if (input_file_number.eq.in_file_number) write(*,*) trim(adjustl(s))

        if (len_trim(s).eq.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        if (index(s, 'f').gt.0 .or. index(s, 'g').gt.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        s = adjustl(remove_tabs(s))
        
        read(s, *, iostat=iostat) array2
        
        if (iostat.eq.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        read(s, *, iostat=iostat) array
        
        if (iostat.ne.0) then
            write(*,*) 'Invalid input, try again.'
            goto 5
        endif
        
        !line_no = line_no + 1
        call write_to_in_file(prompt, array, length)

    end subroutine


    
    subroutine test_prompt(prompt)
		use m_string, only:to_upper
        implicit none
        
        character(*) :: prompt
        character(120) :: temp_string

        integer :: iostat,excl
		logical:: commented
        
        if (input_file_number .ne. 5) then
			commented = .true.
		    do while(commented)
				read(input_file_number, 101, iostat=iostat) temp_string
	101         format(A120)
				!Check for comment character and cycle through file if lines are commented
				excl = index(adjustl(temp_string),'!')
				commented = excl==1
				line_no = line_no + 1
			enddo

			if (excl>0) temp_string = temp_string(:excl-1)

            if (iostat.ne.0) then
                write(*,*) 'Error encountered when reading parameters file.'
                write(*,*) 'End of file reached, but more parameters need to'
                write(*,*) 'be read. Please record the file again.'
                write(*,*)
                pause
                stop
            endif
            
            !line_no = line_no + 1
			!Tabs can cause all sorts of strife so must be removed
            if (to_upper(trim(adjustl(remove_tabs(temp_string)))) .ne. to_upper(trim(prompt))) then
                write(6,*) 'Wrong input string:'
                write(6,*) trim(adjustl(temp_string))
                write(6,*) 'Expected:'
                write(6,*) trim(prompt)
                write(6,100) line_no
100             format(' On line number: ', i3)
                call exit(0)
            endif
        endif
      
    end subroutine



    subroutine init_in_file(fnum)
        
        implicit none
        
        integer(4) :: fnum
        
        output_file_number = fnum
        
    end subroutine

    
    
    subroutine write_prompt(prompt)
        
        implicit none
        
        character(*) :: prompt
        
        write(output_file_number, '(a)') trim(adjustl(prompt))
        
    end subroutine
    
    
    
    subroutine write_to_in_file_int(prompt, num)
    
        implicit none
        
        character(*) :: prompt
        integer(4) :: num
        
        character(1024) :: s
        
        if (output_file_number.gt.0) then
            call write_prompt(prompt)
            
            write(s, *) num
            write(output_file_number, '(4x, a)') trim(adjustl(s))
            
            call flush(output_file_number)
        endif

    end subroutine



    subroutine write_to_in_file_two_ints(prompt, num1, num2)
    
        implicit none
        
        character(*) :: prompt
        integer(4) :: num1, num2
        
        character(1024) :: s
        
        if (output_file_number .gt. 0) then
            call write_prompt(prompt)
            
            write(s, *) num1, num2
            write(output_file_number, '(4x, a)') trim(adjustl(s))
            
            call flush(output_file_number)
        endif

    end subroutine


    subroutine write_to_in_file_int_array(prompt, array, length)
    
        implicit none
        
        character(*) :: prompt
        integer(4) :: length
        integer(4) :: array(length)
        
        character(1024) :: s
        
        if (output_file_number .gt. 0) then
            call write_prompt(prompt)
            
            write(s, *) array
            write(output_file_number, '(4x, a)') trim(adjustl(s))
            
            call flush(output_file_number)
        endif

    end subroutine



    subroutine write_to_in_file_char(prompt, s)
    
        implicit none
        
        character(*) :: prompt
        character(*) :: s
        
        if (output_file_number .gt. 0) then
            call write_prompt(prompt)
            
            write(output_file_number, '(4x, a)') trim(adjustl(s))
            
            call flush(output_file_number)
        endif

    end subroutine


    subroutine write_to_in_file_real(prompt, num)
    
        implicit none
        
        character(*) :: prompt
        real(fp_kind) :: num
        
        character(1024) :: s
        
        if (output_file_number .gt. 0) then
            call write_prompt(prompt)
            
            write(s, *) num
            write(output_file_number, '(4x, a)') trim(adjustl(s))
            
            call flush(output_file_number)
        endif

    end subroutine


    subroutine write_to_in_file_two_reals(prompt, num1, num2)
    
        implicit none
        
        character(*) :: prompt
        real(fp_kind) :: num1, num2
        
        character(1024) :: s
        
        if (output_file_number .gt. 0) then
            call write_prompt(prompt)
            
            write(s, *) num1, num2
            write(output_file_number, '(4x, a)') trim(adjustl(s))
            
            call flush(output_file_number)
        endif

    end subroutine


    subroutine write_to_in_file_real_array(prompt, array, length)
    
        implicit none
        
        character(*) :: prompt
        integer(4) :: length
        real(fp_kind) :: array(length)
        
        character(1024) :: s
        
        if (output_file_number .gt. 0) then
            call write_prompt(prompt)
            
            write(s, *) array
            write(output_file_number, '(4x, a)') trim(adjustl(s))
            
            call flush(output_file_number)
        endif

    end subroutine
	
	function remove_tabs(string)
		use m_string,only:reduce_blanks,replace_text
		implicit none
		character(*),intent(in)::string
		character(:),allocatable :: remove_tabs

		remove_tabs = reduce_blanks(Replace_Text(string,char(9),' '))
	end function

    end module

