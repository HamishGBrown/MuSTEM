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

module m_multislice

    use m_precision, only: fp_kind

    implicit none



    contains

    !This subroutine samples from the available phase grates and then performs one iteration of the multislice algorithm (called in CPU versions only)

    subroutine qep_multislice_iteration(psi,propagator,transmission,nopiy,nopix,ifactory,ifactorx,idum,n_qep_grates,mode&
                                                                    ,shift_arrayy,shift_arrayx,forward_plan,inverse_plan)
        use m_numerical_tools,only:ran1
        use, intrinsic :: iso_c_binding
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in)::propagator(nopiy,nopix),transmission(nopiy,nopix,n_qep_grates)
        complex(fp_kind),intent(in)::shift_arrayy(:,:),shift_arrayx(:,:)
        integer*4,intent(in)::nopiy,nopix,n_qep_grates,mode,ifactory,ifactorx
        type(c_ptr),intent(in),optional::forward_plan,inverse_plan
        integer*4,intent(inout)::idum

        integer*4::shifty,shiftx,nran
        complex(fp_kind)::trans(nopiy,nopix)

        ! Phase grate
        nran = floor(n_qep_grates*ran1(idum)) + 1

        if(mode == 1) then !On the fly calculation
            !call make_qep_potential(trans, tau_slice, nat_slice, ss_slice(7,j))
            !psi_out = psi*trans
        elseif(mode == 2) then !Quick shift
            shiftx = floor(ifactorx*ran1(idum)) * nopix/ifactorx
            shifty = floor(ifactory*ran1(idum)) * nopiy/ifactory
            trans = cshift(cshift(transmission(:,:,nran),shifty,dim=1),shiftx,dim=2)
        elseif(mode == 3) then      !Phase ramp shift
            shiftx = floor(ifactorx*ran1(idum)) + 1
            shifty = floor(ifactory*ran1(idum)) + 1
            call phase_shift_array(transmission(:,:,nran),trans,shift_arrayy(:,shifty),shift_arrayx(:,shiftx))
        else
            trans = transmission(:,:,nran)
        endif
        if(all([present(forward_plan),present(inverse_plan)])) then
            call multislice_iteration(psi,propagator,trans,nopiy,nopix,forward_plan,inverse_plan)
        else
            call multislice_iteration(psi,propagator,trans,nopiy,nopix)
        end if

    end subroutine

    !This subroutine performs one iteration of the multislice algorithm (called in CPU versions only)
    !Probe (psi) input and output is in real space
    subroutine multislice_iteration(psi,propagator,transmission,nopiy,nopix,forward_plan,inverse_plan)
        use FFTW3
        use output
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in)::propagator(nopiy,nopix),transmission(nopiy,nopix)
        integer*4,intent(in)::nopiy,nopix
        type(C_PTR),intent(in),optional :: forward_plan,inverse_plan


        ! Transmit through slice potential (assumes transmission function divided by 1/nopiy/nopix)
        psi = psi*transmission
        ! Propagate to next slice
        if(present(forward_plan)) then
            call inplace_fft(nopiy,nopix,psi,forward_plan)
        else
            call inplace_fft(nopiy,nopix,psi)
        endif
        psi = psi*propagator

        if(present(inverse_plan)) then
            call inplace_ifft(nopiy,nopix,psi,inverse_plan)
        else
            call inplace_ifft(nopiy,nopix,psi)
        endif
    end subroutine

    !This subroutine performs many iterations of the multislice algorithm (called in CPU versions only)
    !Probe (psi) input and output is in real space
    subroutine multislice_iteration_many_slices(psi,propagator,transmission,nopiy,nopix,ncells,nslices,forward_plan,inverse_plan)
        use FFTW3
        use output
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in),dimension(nopiy,nopix,nslices)::propagator,transmission
        integer*4,intent(in)::nopiy,nopix,ncells,nslices
        integer*4::icell,islice
        type(C_PTR),intent(in),optional :: forward_plan,inverse_plan

        do icell=1,ncells;do islice=1,nslices
        ! Transmit through slice potential (assumes transmission function divided by 1/nopiy/nopix)
        psi = psi*transmission(:,:,islice)
        ! Propagate to next slice
        if(present(forward_plan)) then
            call inplace_fft(nopiy,nopix,psi,forward_plan)
        else
            call inplace_fft(nopiy,nopix,psi)
        endif
        psi = psi*propagator(:,:,islice)

        if(present(inverse_plan)) then
            call inplace_ifft(nopiy,nopix,psi,inverse_plan)
        else
            call inplace_ifft(nopiy,nopix,psi)
        endif
        enddo;enddo
    end subroutine
    !This subroutine performs many iterations of the multislice algorithm on the GPU
    !Probe (psi) input and output is in real space
    subroutine multislice_iteration_many_slices_gpu(psi,propagator,transmission,nopiy,nopix,&
                                                    ncells,nslices,forward_plan,inverse_plan)
        use FFTW3
        use output
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in),dimension(nopiy,nopix,nslices)::propagator,transmission
        integer*4,intent(in)::nopiy,nopix,ncells,nslices
        integer*4::icell,islice,y,x
        type(C_PTR),intent(in),optional :: forward_plan,inverse_plan

        !!$acc data
        do icell=1,ncells;do islice=1,nslices
        ! Transmit through slice potential (assumes transmission function divided by 1/nopiy/nopix)
        !!$acc parallel loop
        do y=1,nopiy
            !!$acc parallel loop
            do x =1,nopix
        psi(y,x) = psi(y,x)*transmission(y,x,islice)
        enddo;enddo
        ! Propagate to next slice
        if(present(forward_plan)) then
            call inplace_fft(nopiy,nopix,psi,forward_plan)
        else
            call inplace_fft(nopiy,nopix,psi)
        endif
        !!$acc parallel loop
        do y=1,nopiy
            !!$acc parallel loop
            do x =1,nopix
        psi(y,x) = psi(y,x)*propagator(y,x,islice)
        enddo;enddo

        if(present(inverse_plan)) then
            call inplace_ifft(nopiy,nopix,psi,inverse_plan)
        else
            call inplace_ifft(nopiy,nopix,psi)
        endif
        enddo;enddo
        !!$acc end data
    end subroutine



end module
