      #include "macros.h"
      
      module class_opacity

      use iso_c_binding, only: c_ptr, c_loc, c_f_pointer, c_int, c_double, c_char, c_size_t, c_null_char

      use const_lib
      use crlibm_lib
      use chem_def
      use chem_lib
      use eos_def
      use eos_lib
      use kap_def
      use kap_lib

      implicit none
      private
      public :: Opacity, init_Opacity, shutdown_Opacity, &
            eos_PT, kap_DT, &
            species, num_eos_resuls, &
            nuclide_index

      
      logical, parameter :: use_cache = .true.
      logical, parameter :: use_Type2_opacities = .false.
      integer(c_int), parameter :: species = _SPECIES
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, &
            mg24=7
      character (len=256), parameter :: kappa_file_prefix = 'gn93'
      character (len=256), parameter :: kappa_CO_prefix = 'gn93_co'
      character (len=256), parameter :: kappa_lowT_prefix = &
            'lowT_fa05_gs98'
      
      integer(c_int), protected, bind(C, name="NUM_EOS_RESULTS") :: &
            num_eos_resuls = num_eos_basic_results

      type, bind(C) :: Opacity
         integer(c_int) :: eos_handle, kap_handle
         real(c_double) :: xa(species)
         real(c_double) :: X, Y, Z, abar, zbar, z2bar, ye
         type(c_ptr) :: net_iso, chem_id
      end type Opacity

!      private
!      public :: init_opacity, shutdown_opacity, get_hbar

      contains


      function nuclide_index(c_nuclei) result(indx) bind(C, name='nuclide_index')
         implicit none
         character(c_char), dimension(*), intent(in) :: c_nuclei
         character(len=:), allocatable :: str
         integer(c_int) :: indx
         integer i

         i = 1
         do
            if (c_nuclei(i) == c_null_char) exit
            i = i + 1
         end do

         i = i - 1
         allocate(character(len=i) :: str)
         str = transfer(c_nuclei(1:i), str)

!         type(c_ptr), intent(in) :: c_nuclei
!         character(kind=c_char), pointer, dimension(:) :: nuclei
!         integer(kind=c_size_t), intent(in) :: len_nuclei
!         integer(c_int) :: indx
!
!         call c_f_pointer(c_nuclei, nuclei, [len_nuclei])
         indx = get_nuclide_index(str)
      end function nuclide_index


      subroutine init_const()
         implicit none
         integer :: ierr

         ierr = 0
         call const_init('', ierr)
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            stop 1
         end if
      end subroutine init_const

      subroutine init_chem()
         implicit none
         integer :: ierr

         call crlibm_init
         
         ierr = 0
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in chem_init'
            stop 1
         end if
      end subroutine init_chem


      subroutine init_eos(op)
         implicit none
         type(Opacity), intent(inout) :: op
         integer :: ierr
         real(kind=8) :: frac, dabar_dx(species), dzbar_dx(species), &
               sumx, xh, xhe, xz, mass_correction, dmc_dx(species)
         integer(c_int), pointer, dimension(:) :: net_iso, chem_id

         call eos_init('mesa', '', '', '', use_cache, ierr)
         if (ierr /= 0) then
            write(*,*) 'eos_init failed'
            stop 1
         end if
         
         op%eos_handle = alloc_eos_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed trying to allocate eos eos_handle'
            stop 1
         end if
         
         allocate(net_iso(num_chem_isos), chem_id(species), &
               stat=ierr)
         op%net_iso = c_loc(net_iso(1))
         op%chem_id = c_loc(chem_id(1))
         if (ierr /= 0) stop 'allocate failed'
         net_iso(:) = 0
         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         op%X = op%xa(h1)
         op%Y = op%xa(he4)
         op%Z = 1 - (op%X + op%Y)
         call composition_info( &
               species, chem_id, op%xa, xh, xhe, xz, op%abar, &
               op%zbar, op%z2bar, op%ye, mass_correction, sumx, &
               dabar_dx, dzbar_dx, dmc_dx)
      end subroutine init_eos


      subroutine init_kap(op)
         implicit none
         type(Opacity), intent(inout) :: op
         integer :: ierr

         call kap_init(kappa_file_prefix, kappa_CO_prefix, &
               kappa_lowT_prefix, 0.0_dp, 0.0_dp, use_cache, &
               '', '', ierr)
         if(ierr/=0) stop 'problem in kap_init'
         
         op%kap_handle = alloc_kap_handle(ierr)
         if(ierr/=0) stop 'problem in alloc_kap_handle'

         ! All these numbers do not matter while use_Type2_opacity is false
         call kap_set_choices(op%kap_handle, .true., .true., .true., &
               .true., use_Type2_opacities, &
               0.71_dp, 0.70_dp, 0.001_dp, 0.01_dp, &
               ierr)
         if(ierr/=0) stop 'problem in kap_set_interpolation_choices'
      end subroutine init_kap


      subroutine init_Opacity(op) bind(C, name='init_Opacity')
         implicit none
         type(Opacity), intent(inout) :: op
         
         call init_const
         call init_chem
         call init_eos(op)
         call init_kap(op)
      end subroutine init_opacity


      subroutine shutdown_eos(op)
         implicit none      
         type(Opacity), intent(inout) :: op
         !integer, pointer, dimension(:) :: net_iso, chem_id

         call free_eos_handle(op%eos_handle)
         call eos_shutdown
         ! call c_f_pointer(op%net_iso, net_iso, [num_chem_isos])
         ! call c_f_pointer(op%chem_id, chem_id, [species])
         ! deallocate(net_iso, chem_id)
      end subroutine shutdown_eos


      subroutine shutdown_kap(op)
         implicit none
         type(Opacity), intent(inout) :: op

         call free_kap_handle(op%kap_handle)
         call kap_shutdown
      end subroutine shutdown_kap


      subroutine shutdown_Opacity(op) bind(C, name='shutdown_Opacity')
         type(Opacity) :: op
         
         call shutdown_eos(op)
         call shutdown_kap(op)
      end subroutine shutdown_opacity

      
      subroutine eos_PT(op, Pgas, T, &
            Rho, log10Rho, dlnRho_dlnPgas_const_T, &
            dlnRho_dlnT_const_Pgas, res, &
            ierr &
            ) bind(C, name='eos_PT')
         implicit none
         type(Opacity), intent(in) :: op
         real(c_double), value :: Pgas, T
         real(c_double), intent(out) :: Rho, log10Rho, &
               dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
         real(c_double), intent(inout) :: res(num_eos_basic_results)
         integer(c_int), intent(out) :: ierr
         real(c_double), dimension(num_eos_basic_results) :: &
               d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_const_TRho, d_dzbar_const_TRho
         integer(c_int), pointer, dimension(:) :: net_iso, chem_id

         call c_f_pointer(op%net_iso, net_iso, [num_chem_isos])
         call c_f_pointer(op%chem_id, chem_id, [species])
         call eosPT_get(op%eos_handle, op%Z, op%X, op%abar, op%zbar, &
               species, chem_id, net_iso, op%xa, &
               Pgas, log10_cr(Pgas), T, log10_cr(T), &
               Rho, log10Rho, &
               dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_const_TRho, d_dzbar_const_TRho, &
               ierr)
      end subroutine eos_PT


      subroutine kap_DT(op, Rho, T, lnfree_e, &
            kappa, dlnkap_dlnRho, dlnkap_dlnT, ierr &
            ) bind(C, name='kap_DT')
         implicit none
         type(Opacity), intent(in) :: op
         real(c_double), value, intent(in) :: Rho, T, lnfree_e
         real(c_double), intent(out) :: kappa, dlnkap_dlnRho, dlnkap_dlnT
         integer(c_int), intent(out) :: ierr
         real(kind=8), parameter :: d_lnfree_e_dlnRho = 0.0  ! needed for Compton
         real(kind=8), parameter :: d_lnfree_e_dlnT = 0.0  ! needed for Compton
         real(kind=8) :: frac_Type2

         

         call kap_get(op%kap_handle, &
               op%zbar, op%X, op%Z, op%Z, &
               op%xa(c12), op%xa(n14), op%xa(o16), op%xa(ne20), &
               log10_cr(Rho), log10_cr(T), &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               frac_Type2, kappa, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      end subroutine kap_DT


      end module class_opacity

