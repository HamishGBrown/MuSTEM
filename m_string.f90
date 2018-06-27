module m_string
use m_precision
implicit none
	interface to_string
		module procedure to_string_integer,to_string_real
	end interface
    

    contains
      
	function logical_to_yn(input) result(str)
		logical,intent(in)::input
		character(1)::str
		if (input) str = 'y'
		if (.not.input) str = 'n'
    end function
    
		 function is_numeric(string)
			  character(len=*), intent(in) :: string
			  logical :: is_numeric
			  real :: x
			  integer :: e
			  read(string,*,iostat=e) x
			  is_numeric = e == 0
			end function is_numeric
			
        function to_string_integer(i) result(s)
            implicit none

            character(:),allocatable :: s

            integer :: i
            
            character(32) :: ss
            integer :: l
            
            write(ss, *) i
            ss = adjustl(ss)

            l = len_trim(ss)

            allocate(character(l)::s)

            s = trim(ss)

        end function
			
    pure function to_string_real(input) result(to_string)
        real(fp_kind),intent(in):: input
        character(len=:),allocatable :: to_string
        character*50::dum
        if((abs(input).lt.1d3).and.(abs(input).gt.1d-3)) then
			write(dum, '(f9.3)') input
		elseif((abs(input).lt.1d-16)) then
			write(dum, '(i1)') 0
		else
			write(dum, '(e10.3)') input
		endif
        allocate(character(len = len_trim(adjustl(dum))) :: to_string)
        to_string = trim(adjustl(dum))
    end function
    
    elemental function str2int(str) result(int)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer         :: int
    integer :: stat_

    read(str,*,iostat=stat_)  int
  end function str2int
    
	!Allows thickness to be read in via:
	!Single value in  the format: [z1]
	!List of values in the format: [z1],[z2],[z3]
	!A sequence of values in the format: [zstart]:[zstop]:[zstep]
	subroutine read_sequence_string(string,lstring,nz,zarray,minstep)
		character(len=lstring), intent(in) :: string
		integer*4,intent(in)::lstring
		integer*4,intent(inout)::nz
		real(fp_kind),optional::zarray(nz),minstep
		
		character(len=1)::sgn
		integer*4::i,ii,j,nsteps,zi
		real(fp_kind):: zseq(3),z
		logical:: is_sequence,is_list
		
		
		!See if ':' is string
		is_sequence = (index(string,':').ne.0)
		
		!See if ',' is string
		is_list = (index(string,',').ne.0)
		
		if (is_sequence.and.is_list) then
			pause "Both of the characters ',', for a list of values, and ':', for an ordered sequence of values, were found in the input, please rectify and re-run program"
			stop
		endif

		if(is_sequence) then
			i=1
			do ii=1,2
				j = index(string(i:),':')+i-2

				if (is_numeric(string(i:j))) then
					read(string(i:j),*) zseq(ii)
				else
					write(*,*) "Expected numeric input of the kind [start]:[stop]:[step], got: "//string(i:j)
					stop
				endif
				i=j+2
			enddo

			!Read final value into array
			if (is_numeric(string(i:))) then
				read(string(i:),*) zseq(ii)
				if(present(minstep)) zseq(ii) = max(minstep,zseq(ii))
			else
				write(*,*) "Expected numeric input of the kind [start]:[stop]:[step], got: "//string(i:)
				stop
			endif

			if(zseq(2)<zseq(1)*sign(1.0_fp_kind,zseq(3))) then
                if(zseq(3)>0) then
                    sgn = ">"
                else
                    sgn = "<"
                endif
				write(*,*) "Expected numeric input of the kind [start]:[stop]:[step], however [start] = "//to_string(zseq(1))//sgn//" [stop] = "//to_string(zseq(2))//" with [step] = "//to_string(zseq(3))
				stop
			endif
			z=zseq(1)
			nsteps = 0
			do while(z<=zseq(2))
				nsteps = nsteps+1
				z = z+zseq(3)
			enddo
			!nsteps = nint((zseq(2)-zseq(1)/zseq(3)))+1

			if (present(zarray)) then
				zarray = (/(zi, zi=0,nsteps, 1)/)*zseq(3)+zseq(1)
			else
				nz = nsteps
			endif
			
		elseif(is_list) then 
			i=1
			ii=0
			do while (index(string(i:),',').ne.0)  
				!write(*,*) i,index(string(i:),',')
				j = index(string(i:),',')+i-2
				ii=ii+1
				if (present(zarray)) then 
					if (.not.is_numeric(string(i:j))) then
						write(*,*) "Expected numeric input, got: "//string(i:j)
						stop
					endif
					read(string(i:j),*) zarray(ii)
				endif
				
				i=j+2
				
				
			end do
			nz = ii+1
			!Read final value into array
			if (.not.is_numeric(string(i:))) then
				write(*,*) "Expected numeric input, got: "//string(i:)
				stop
			endif
			if (present(zarray)) read(string(i:),*) zarray(ii+1)
		else
			
			if (present(zarray)) then
				
				read(string,*) zarray(1)
			else
				
				nz = 1
				
			endif
		endif
	
	end subroutine

	function zero_padded_int(num, width)
        implicit none
        
        integer*4,intent(in) :: num, width
        character(len=width) :: zero_padded_int
        
        
        character(10) :: fmt
        
10      format('(i', i1, '.', i1, ')')   

        if (num.ge.0) then
            write(fmt, 10) width, width
        else
            write(fmt, 10) width, width-1
        endif

        write(zero_padded_int, fmt) num
        
      end function

	function to_upper(strIn) result(strOut)
		! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
		! Original author: Clive Page

			 implicit none

			 character(len=*), intent(in) :: strIn
			 character(len=len(strIn)) :: strOut
			 integer :: i,j

			 do i = 1, len(strIn)
				  j = iachar(strIn(i:i))
				  if (j>= iachar("a") .and. j<=iachar("z") ) then
					   strOut(i:i) = achar(iachar(strIn(i:i))-32)
				  else
					   strOut(i:i) = strIn(i:i)
				  end if
			 end do

	end function to_upper

	function to_lower(strIn) result(strOut)
		! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
		! Original author: Clive Page

			 implicit none

			 character(len=*), intent(in) :: strIn
			 character(len=len(strIn)) :: strOut
			 integer :: i,j

			 do i = 1, len(strIn)
				  j = iachar(strIn(i:i))
				  if (j>= iachar("A") .and. j<=iachar("Z") ) then
					   strOut(i:i) = achar(iachar(strIn(i:i))+32)
				  else
					   strOut(i:i) = strIn(i:i)
				  end if
			 end do

	end function to_lower

	FUNCTION Reduce_Blanks(s)  RESULT (outs)
		CHARACTER(*)      :: s
		CHARACTER(LEN_TRIM(s)) :: outs
		INTEGER           :: i, k, n

		n = 0  ; k = LEN_TRIM(s)          ! k=index last non-blank (may be null)
		DO i = 1,k-1                      ! dont process last char yet
		   n = n+1 ; outs(n:n) = s(i:i)
		   IF (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
		END DO
		n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
		IF (n < k) outs(n+1:) = ' '       ! pad trailing blanks
		END FUNCTION Reduce_Blanks

	! ------------------
	FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
		CHARACTER(*)        :: s,text,rep
		CHARACTER(LEN(s)+100) :: outs     ! provide outs with extra 100 char len
		INTEGER             :: i, nt, nr

		outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
		DO
		   i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
		   outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
		END DO
	END FUNCTION Replace_Text
end module

