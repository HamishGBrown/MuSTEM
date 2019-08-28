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

    module global_variables
    
    use m_precision, only: fp_kind
    
    implicit none
    
    integer(4) :: nt,nm,i_xtl                   !number of atom types, number of atoms, xtl file read flag
    integer(4) :: nopiy,nopix,npixels           !supercell
    integer(4) :: nopix_ucell,nopiy_ucell       !unit cell
    integer(4) :: ifactorx,ifactory             !unit cell tilings
    real(fp_kind) :: deltay,deltax              !real space spacing between pixels in x and y
    real(fp_kind) :: normalisation
    
    real(fp_kind), allocatable :: bwl_mat(:,:)            !bandwidth limiting matrix
                                          
    integer(4), allocatable    :: nat(:)                     !number of each atom type in the unit celll
	real(fp_kind), allocatable	:: dz(:)					  !ionicity of each atom type
    real(fp_kind), allocatable :: tau(:,:,:)              !position of the atoms in the unit cell
    real(fp_kind), allocatable :: atf(:,:)                !atomic number, occupancy and DWF (urms)
    real(fp_kind), allocatable :: atomf(:,:),fx(:)        !electron scattering factor parameterisation from elsa
    real(fp_kind)  :: a0(3),deg(3),ekv,ss(7)              !a b c unit cell lengths, angle between a b c, accelerating voltage, tricyclinc info
    real(fp_kind)  :: thetad,surfn(3),orthog(3,3)
    real(fp_kind)  :: volts,ak                            !mean inner potential, wavevector (corrected for refraction)
    real(fp_kind)  :: ak1,relm                            !wavevector in freespace, relativistically corrected mass
    real(fp_kind),allocatable  :: claue(:,:),Kz(:)        !Specimen tilt vector and z component of incident wave vector
    integer*4::n_tilts_total                              !Total number of specimen tilts
    real(fp_kind)  :: bvec(3)                             !Beam tilt vector
    
    !sample thickness and slicing variables
    real(fp_kind) :: thickness                        !sample thickness
	real(fp_kind),allocatable:: zarray(:)
	integer(4),allocatable::ncells(:)
	integer(4)::nz
    integer(4) :: n_cells  

    complex(fp_kind), allocatable :: fz(:,:,:)            !the scattering factors, in reciprocal space, calculated on the grid (supercell)
    complex(fp_kind), allocatable :: fz_DWF(:,:,:)        !the DWF smear_array, in reciprocal space, calculated on the grid (supercell)
    complex(fp_kind), allocatable :: sinc(:,:)            !sinc function to correct for pixelation in the potential construction
    complex(fp_kind), allocatable :: inverse_sinc(:,:)    !1/sinc function to correct for pixelation in the potential construction

    real(fp_kind)  :: uvw1(3),uvw2(3)                     !real space scan vectors that are parallel 
    integer(4) :: ig1(3),ig2(3),izone(3)                  !to the reciprocal space vectors ig1,ig2.
        
    character*120 :: substance                            !label for the crystal substance
    character*10, allocatable :: substance_atom_types(:)

    !output variables
    integer(4) :: ndet,nseg ,nopiyout,nopixout                           !number of integrating detectors and 4D STEM output
    real(fp_kind) :: seg_det_offset
	logical::segments
    real(fp_kind), allocatable :: outer(:),inner(:)       !detector ranges (in inverse angstrom)
    
    !interpolation variables
    integer(4) :: output_nopiy,output_nopix               !output number of pixels in the interpolated output image
    integer(4) :: tilex,tiley                             !Interpolated image tiling
    real(fp_kind)  ::  bwl_rad                            !band width limiting radius (in q space)
    integer(4) :: bin_factor                              !Factor by which to bin 4DSTEM diffraction pattern outputs
                                                
    !Constants data
    real(fp_kind),parameter :: pi = atan(1.0_fp_kind)*4.0_fp_kind 
    real(fp_kind),parameter :: tp = atan(1.0_fp_kind)*8.0_fp_kind 
    real(fp_kind),parameter :: hsq_on_twom = 150.4132_fp_kind  !h^2/2m
    complex(fp_kind) :: ci = cmplx(0.0_fp_kind,1.0_fp_kind)
    
    real(fp_kind),parameter :: const1 = 9.7846113e-07_fp_kind       ! const1 = 1/(2mc^2) in eV^(-1)
    real(fp_kind),parameter :: const2 = 12.263868_fp_kind           ! const2 = h  in eV.A
    real(fp_kind),parameter :: bohr = 0.529177_fp_kind              ! bohr radius
    real(fp_kind),parameter :: ryd = 13.60535_fp_kind               ! rydberg constant
    real(fp_kind),parameter :: fsc = 7.29735e-3_fp_kind             ! fsc = fine structure constant (dimensionless)
    real(fp_kind),parameter :: hbarc = 1973.26_fp_kind              ! hbarc = hbar * c in eV A units
    
    !logical types to pick inelastic calculations
    logical :: adf 
    logical :: EELS = .false.

    
    logical :: qep,output_thermal,interpolation,fourdSTEM
    
    logical :: on_the_fly = .false.
    logical :: high_accuracy
	logical :: ionic = .false.
    
	contains
     
    
      function wavev(e)
      !  this function returns the wavevector in one over lambda, in a-1,
      !  for an input electron energy e in ev.
      
      use m_precision
      
      implicit none
      
      real(fp_kind) c1,c2,e
      real(fp_kind) wavev
      data c1, c2 / 9.78475598e-07_fp_kind, 12.2642596_fp_kind /
      
      wavev = sqrt( e + c1 *e ** 2.0_fp_kind ) / c2
      
      end function 
    
      subroutine constants()
      
          use m_precision
      
          implicit none
      
      
          real(fp_kind) con1,c2,delk
          data con1, c2 / 510.9989461_fp_kind, 1.956951198e-03_fp_kind /

          relm = ( con1 + ekv) * c2
          
    ! ak is the incident wave vector in the solid corrected for refraction (mean inner potential)
          ak   = wavev( ekv * 1000_fp_kind + volts)
          
    ! ak1 is the incident wave vector without corrections for refraction
          ak1   = wavev( ekv * 1000_fp_kind )
    !Initialise tilt to zero (on-axis)
          allocate(Kz(1),claue(3,1))
          n_tilts_total = 1
          Kz = ak1
          claue = 0_fp_kind
          
          delk = ak - ak1

	    write(6,*) 'Pertinent quantities for the crystal:'
          write(6,111) ekv,ak,char(143),ak1,char(143),delk,char(143),relm
      111 format('          E = ',F12.3,' keV',/,&
         &' Crystal Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'  Vacuum Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'   Delta Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'     m* / m = ',g12.5,&
         &' (relativistic mass increment)',/,/)

      end subroutine 
      

      end module global_variables
