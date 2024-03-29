
      module class_opacity

      use iso_c_binding, only: c_ptr, c_loc, c_f_pointer, c_int, c_double, c_char, c_size_t, c_null_char

      use const_lib
      use chem_def
      use chem_lib
      use eos_def
      use eos_lib
      use kap_def
      use kap_lib
      use math_lib

      implicit none
      private
      public :: init_mesa, get_num_chem_isos, &
            Opacity, init_Opacity, shutdown_Opacity, &
            eos_PT, kap_DT, num_eos_resuls, &
            nuclide_index, solsiz, solx, get_sol_x

      logical, parameter :: use_cache = .true.
      character (len=256), parameter :: kappa_file_prefix = 'gs98'
      character (len=256), parameter :: kappa_CO_prefix = 'gs98_co'
      character (len=256), parameter :: kappa_lowT_prefix = &
            'lowT_fa05_gs98'
      
      integer(c_int), protected, bind(C, name="NUM_EOS_RESULTS") :: &
            num_eos_resuls = num_eos_basic_results

      integer(c_int), protected, bind(C, name="SOLSIZE") :: &
            solsize = solsiz

      type, bind(C) :: Opacity
         integer(c_int) :: eos_handle, kap_handle
         integer(c_int) :: species
         real(c_double) :: X, Y, Z, XC, XN, XO, XNe, abar, zbar, z2bar, ye
         type(c_ptr) :: net_iso, chem_id, xa
      end type Opacity

!      private
!      public :: init_opacity, shutdown_opacity, get_hbar

      contains


      subroutine get_sol_x(c_sol_x) bind(C, name='get_sol_x')
          type(c_ptr), value :: c_sol_x
          real(c_double), pointer, dimension(:) :: sol_x
          call c_f_pointer(c_sol_x, sol_x, [solsiz])
          sol_x = solx
      end subroutine


      subroutine get_sol_chem_id(c_chem_id) bind(C, name='get_sol_chem_id')
          integer i
          type(c_ptr), value :: c_chem_id
          integer(c_int), pointer, dimension(:) :: chem_id
          call c_f_pointer(c_chem_id, chem_id, [solsiz])
          do i=1, solsiz
              chem_id(i) = get_nuclide_index(namsol(i))
          end do
      end subroutine get_sol_chem_id


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


      subroutine mesa_init_const()
         implicit none
         integer :: ierr

         ierr = 0
         call const_init('', ierr)
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            stop 1
         end if
      end subroutine mesa_init_const

      subroutine mesa_init_math()
         call math_init()
      end subroutine mesa_init_math

      subroutine mesa_init_chem()
         implicit none
         integer :: ierr

         ierr = 0
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'chem_init failed'
            stop 1
         end if
      end subroutine mesa_init_chem

      subroutine mesa_init_eos()
         implicit none
         integer :: ierr
         
         call eos_init('mesa', '', '', '', use_cache, ierr)
         if (ierr /= 0) then
            write(*,*) 'eos_init failed'
            stop 1
         end if
      end subroutine mesa_init_eos

      subroutine mesa_init_kap()
         implicit none
         integer :: ierr
         real(kind=8), parameter :: kap_blend_logT_upper_bdy = 3.88_dp  ! default
         real(kind=8), parameter :: kap_blend_logT_lower_bdy = 3.80_dp  ! default

         call kap_init(kappa_file_prefix, kappa_CO_prefix, &
                 kappa_lowT_prefix, kap_blend_logT_upper_bdy, &
                 kap_blend_logT_lower_bdy, use_cache, &
                 '', '', .false., ierr)
         if(ierr /= 0) then
            write(*,*) 'kap_init failed'
            stop 1
         end if
      end subroutine mesa_init_kap

      subroutine init_mesa() bind(C, name='init_mesa')
         call mesa_init_const()
         call mesa_init_math()
         call mesa_init_chem()
         call mesa_init_eos()
         call mesa_init_kap()
      end subroutine init_mesa

      subroutine mesa_shutdown_eos()
         call eos_shutdown()
      end subroutine mesa_shutdown_eos

      subroutine mesa_shutdown_kap()
         call kap_shutdown()
      end subroutine mesa_shutdown_kap

      subroutine shutdown_mesa() bind(C, name='shutdown_mesa')
         call mesa_shutdown_eos()
         call mesa_shutdown_kap()
      end subroutine shutdown_mesa


      subroutine init_eos(op)
         implicit none
         type(Opacity), intent(inout) :: op
         integer :: ierr, i
         real(kind=8) :: frac, dabar_dx(op%species), dzbar_dx(op%species), &
               sumx, xh, xhe, xz, mass_correction, dmc_dx(op%species)
         integer(c_int), pointer, dimension(:) :: chem_id, net_iso
         real(c_double), pointer, dimension(:) :: xa
         integer :: ihe6, ihe7, ihe8, ihe9, ihe10
         integer, dimension(3) :: ih_array
         integer, dimension(7) :: ihe_array
         integer, dimension(8) :: ic_array
         integer, dimension(9) :: in_array
         integer, dimension(8) :: io_array
         integer, dimension(12) :: ine_array

         call c_f_pointer(op%chem_id, chem_id, [op%species])
         call c_f_pointer(op%net_iso, net_iso, [num_chem_isos])
         call c_f_pointer(op%xa, xa, [op%species])
         
         op%eos_handle = alloc_eos_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed trying to allocate eos eos_handle'
            stop 1
         end if

         op%X = 0
         op%Y = 0
         op%XC = 0
         op%XN = 0
         op%XO = 0
         op%XNe = 0

         ihe6 = get_nuclide_index('he6')
         ihe7 = get_nuclide_index('he7')
         ihe8 = get_nuclide_index('he8')
         ihe9 = get_nuclide_index('he9')
         ihe10 = get_nuclide_index('he10')

         ih_array = (/ ih1, ih2, ih3 /)
         ihe_array = (/ ihe3, ihe4, ihe6, ihe7, ihe8, ihe9, ihe10 /)
         ic_array = (/ ic9, ic10, ic11, ic12, ic13, ic14, ic15, ic16 /)
         in_array = (/ in12, in13, in14, in15, in16, in17, in18, in19, in20 /)
         io_array = (/ io13, io14, io15, io16, io17, io18, io19, io20 /)
         ine_array = (/ ine17, ine18, ine19, ine20, ine21, ine22, ine23, ine24, ine25, ine26, ine27, ine28 /)

         do i = 1, size(ih_array)
             if (net_iso(ih_array(i)) /= 0) then
                 op%X = op%X + xa(net_iso(ih_array(i)))
             end if
         end do

         do i = 1, size(ihe_array)
             if (net_iso(ihe_array(i)) /= 0) then
                 op%Y = op%Y + xa(net_iso(ihe_array(i)))
             end if
         end do

         do i = 1, size(ic_array)
             if (net_iso(ic_array(i)) /= 0) then
                 op%XC = op%XC + xa(net_iso(ic_array(i)))
             end if
         end do

         do i = 1, size(in_array)
             if (net_iso(in_array(i)) /= 0) then
                 op%XN = op%XN + xa(net_iso(in_array(i)))
             end if
         end do

         do i = 1, size(io_array)
             if (net_iso(io_array(i)) /= 0) then
                 op%XO = op%XO + xa(net_iso(io_array(i)))
             end if
         end do

         do i = 1, size(ine_array)
             if (net_iso(ine_array(i)) /= 0) then
                 op%XNe = op%XNe + xa(net_iso(ine_array(i)))
             end if
         end do

         op%Z = 1 - (op%X + op%Y)

         if (op%Z /= 0) then
             op%XC = op%XC / op%Z
             op%XN = op%XN / op%Z
             op%XO = op%XO / op%Z
             op%XNe = op%XNe / op%Z
         end if

         call composition_info( &
               op%species, chem_id, xa, xh, xhe, xz, op%abar, &
               op%zbar, op%z2bar, op%ye, mass_correction, sumx, &
               dabar_dx, dzbar_dx, dmc_dx)
      end subroutine init_eos


      subroutine init_kap(op)
         implicit none
         type(Opacity), intent(inout) :: op
         integer :: ierr

         logical, parameter :: cubic_interpolation_in_X = .true.
         logical, parameter :: cubic_interpolation_in_Z = .true.
         logical, parameter :: include_electron_conduction = .true.
         ! We are not going to use Zbase which is for CN abundance
         logical, parameter :: use_Zbase_for_Type1 = .false.
         ! We don not need Type2 opacities, which are useful
         ! for CN-rich composition only.
         ! All Type2 variables are not used
         logical, parameter :: use_Type2_opacities = .false.
         real(kind=8), parameter :: kap_Type2_full_off_X = 0.71_dp
         real(kind=8), parameter :: kap_Type2_full_on_X = 0.70_dp
         real(kind=8), parameter :: kap_Type2_full_off_dZ = 0.001_dp
         real(kind=8), parameter :: kap_Type2_full_on_dZ = 0.01_dp


         op%kap_handle = alloc_kap_handle(ierr)
         if(ierr/=0) stop 'problem in alloc_kap_handle'


         call kap_set_choices( &
               op%kap_handle, cubic_interpolation_in_X, cubic_interpolation_in_Z, &
               include_electron_conduction, &
               use_Zbase_for_Type1, use_Type2_opacities, &
               kap_Type2_full_off_X, kap_Type2_full_on_X, &
               kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
               ierr)

         if(ierr/=0) stop 'problem in kap_set_interpolation_choices'
      end subroutine init_kap

      function get_num_chem_isos() result(n) bind(C, name='get_num_chem_isos')
          implicit none
          integer(c_int) :: n

          n = num_chem_isos
      end function get_num_chem_isos

      subroutine init_Opacity(op) bind(C, name='init_Opacity')
         implicit none
         type(Opacity), intent(inout) :: op

         call init_eos(op)
         call init_kap(op)
      end subroutine init_opacity


      subroutine shutdown_eos(op)
         implicit none      
         type(Opacity), intent(inout) :: op
         ! integer, pointer, dimension(:) :: net_iso, chem_id

         call free_eos_handle(op%eos_handle)
         ! call c_f_pointer(op%net_iso, net_iso, [num_chem_isos])
         ! call c_f_pointer(op%chem_id, chem_id, [species])
         ! deallocate(net_iso, chem_id)
      end subroutine shutdown_eos


      subroutine shutdown_kap(op)
         implicit none
         type(Opacity), intent(inout) :: op

         call free_kap_handle(op%kap_handle)
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
         real(c_double), pointer, dimension(:) :: xa

         call c_f_pointer(op%net_iso, net_iso, [num_chem_isos])
         call c_f_pointer(op%chem_id, chem_id, [op%species])
         call c_f_pointer(op%xa, xa, [op%species])


         call eosPT_get(op%eos_handle, op%Z, op%X, op%abar, op%zbar, &
               op%species, chem_id, net_iso, xa, &
               Pgas, log10(Pgas), T, log10(T), &
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

         ! We do not need Type2 opacities which use Zbase for CN abundance
         real(kind=8), parameter :: Zbase = -1.0_dp

         call kap_get(op%kap_handle, &
               op%zbar, op%X, op%Z, Zbase, &
               op%XC, op%XN, op%XO, op%XNe, &
               log10(Rho), log10(T), &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               frac_Type2, kappa, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      end subroutine kap_DT


      end module class_opacity
