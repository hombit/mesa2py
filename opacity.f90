      module class_opacity

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
            eos_PT, kap_log10DT, &
            species

      
      logical, parameter :: use_cache = .true.
      logical, parameter :: use_Type2_opacities = .false.
      integer, parameter :: species = 7
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, &
            mg24=7
      
      
      type :: Opacity
         integer :: eos_handle, kap_handle
         real(kind=8) :: xa(species)
         real(kind=8) :: Y, abar, zbar, z2bar, ye
         integer, pointer, dimension(:) :: net_iso, chem_id
         real(kind=8) :: X = 0.70
         real(kind=8) :: Z = 0.02
         real(kind=8) :: Zfrac_C = 0.173312d0
         real(kind=8) :: Zfrac_N = 0.053177d0
         real(kind=8) :: Zfrac_O = 0.482398d0
         real(kind=8) :: Zfrac_Ne = 0.098675d0
      end type Opacity

!      private
!      public :: init_opacity, shutdown_opacity, get_hbar

      contains


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
         character (len=256) :: eos_file_prefix
         integer :: ierr
         integer, pointer, dimension(:) :: net_iso, chem_id
         real(kind=8) :: frac, dabar_dx(species), dzbar_dx(species), &
               sumx, xh, xhe, xz, mass_correction, dmc_dx(species)

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
         
         allocate(op%net_iso(num_chem_isos), op%chem_id(species), &
               stat=ierr)
         if (ierr /= 0) stop 'allocate failed'
         op%net_iso(:) = 0
         op%chem_id(h1) = ih1; op%net_iso(ih1) = h1
         op%chem_id(he4) = ihe4; op%net_iso(ihe4) = he4
         op%chem_id(c12) = ic12; op%net_iso(ic12) = c12
         op%chem_id(n14) = in14; op%net_iso(in14) = n14
         op%chem_id(o16) = io16; op%net_iso(io16) = o16
         op%chem_id(ne20) = ine20; op%net_iso(ine20) = ne20
         op%chem_id(mg24) = img24; op%net_iso(img24) = mg24
         op%Y = 1 - (op%X + op%Z)               
         op%xa(h1) = op%X
         op%xa(he4) = op%Y
         op%xa(c12) = op%Z * op%Zfrac_C
         op%xa(n14) = op%Z * op%Zfrac_N
         op%xa(o16) = op%Z * op%Zfrac_O
         op%xa(ne20) = op%Z * op%Zfrac_Ne
         op%xa(species) = 1 - sum(op%xa(1:species-1))
         call composition_info( &
               species, op%chem_id, op%xa, xh, xhe, xz, op%abar, &
               op%zbar, op%z2bar, op%ye, mass_correction, sumx, &
               dabar_dx, dzbar_dx, dmc_dx)
      end subroutine init_eos


      subroutine init_kap(op)
         implicit none
         type(Opacity), intent(inout) :: op
         integer :: ierr

         call kap_init('gn93', '', '', 0.0_dp, 0.0_dp, use_cache, &
               '', '', ierr)
         if(ierr/=0) stop 'problem in kap_init'
         
         op%kap_handle = alloc_kap_handle(ierr)
         if(ierr/=0) stop 'problem in alloc_kap_handle'

         ! All these numbers don't matter while use_Type2_opacity is false
         call kap_set_choices(op%kap_handle, .true., .true., .true., &
               .true., use_Type2_opacities, &
               0.71_dp, 0.70_dp, 0.001_dp, 0.01_dp, &
               ierr)
         if(ierr/=0) stop 'problem in kap_set_interpolation_choices'
      end subroutine init_kap


      type(Opacity) function init_opacity()
         implicit none
         
         call init_const
         call init_chem
         call init_eos(init_opacity)
         call init_kap(init_opacity)
      end function init_opacity


      subroutine shutdown_eos(op)
         implicit none      
         type(Opacity), intent(inout) :: op

         call free_eos_handle(op%eos_handle)
         call eos_shutdown
         deallocate(op%net_iso, op%chem_id)
      end subroutine shutdown_eos


      subroutine shutdown_kap(op)
         implicit none
         type(Opacity), intent(inout) :: op

         call free_kap_handle(op%kap_handle)
         call kap_shutdown
      end subroutine shutdown_kap


      subroutine shutdown_opacity(op)
         type(Opacity) :: op
         
         call shutdown_eos(op)
         call shutdown_kap(op)
      end subroutine shutdown_opacity
           
      
      subroutine eos_PT(op, Pgas, T, &
            Rho, log10Rho, dlnRho_dlnPgas_const_T, &
            dlnRho_dlnT_const_Pgas, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_const_TRho, d_dzbar_const_TRho, &
            ierr)
         implicit none
         type(Opacity), intent(in) :: op
         real(kind=8), intent(in) :: Pgas, T
         real(kind=8), intent(out) :: Rho, log10Rho, &
               dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
         real(kind=8), intent(inout) :: res(:)
         real(kind=8), intent(inout) :: d_dlnRho_const_T(:), &
               d_dlnT_const_Rho(:)
         real(kind=8), intent(inout) :: d_dabar_const_TRho(:), &
               d_dzbar_const_TRho(:)
         integer, intent(out) :: ierr

         call eosPT_get(op%eos_handle, op%Z, op%X, op%abar, op%zbar, &
               species, op%chem_id, op%net_iso, op%xa, &
               Pgas, log10_cr(Pgas), T, log10_cr(T), &
               Rho, log10Rho, &
               dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_const_TRho, d_dzbar_const_TRho, &
               ierr)
      end subroutine eos_PT


      subroutine kap_log10DT(op, log10Rho, log10T, &
            kappa, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         implicit none
         type(Opacity), intent(in) :: op
         real(kind=8), intent(in) :: log10Rho, log10T
         real(kind=8), intent(out) :: kappa, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr
         real(kind=8), parameter :: lnfree_e = 0.0  ! needed for Compton
         real(kind=8), parameter :: d_lnfree_e_dlnRho = 0.0
         real(kind=8), parameter :: d_lnfree_e_dlnT = 0.0
         real(kind=8) :: frac_Type2

         call kap_get(op%kap_handle, &
               op%zbar, op%X, op%Z, op%Z, &
               op%xa(c12), op%xa(n14), op%xa(o16), op%xa(ne20), &
               log10Rho, log10T, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               frac_Type2, kappa, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      end subroutine kap_log10DT


      end module class_opacity

