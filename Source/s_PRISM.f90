subroutine multislice_calc(nopiy,nopix,ifactory,ifactorx)

    complex(fp_kind)::psi_initial(nopiy,nopix)

    !Set up illumination
    if(pw_illum)
        psi_initial = 1.0_fp_kind/sqrt(float(nopiy*nopix))
    else
        psi_initial = make_ctf(probe_initial_position,probe_df(1),probe_cutoff,probe_aberrations,probe_apodisation)
        call inplace_ifft(nopiy, nopix, psi_initial)
        psi_initial = psi_initial/sqrt(sum(abs(psi_initial)**2))
    endif


    !Loop over frozen phonon configurations
    do ifp=1,
end subroutine
