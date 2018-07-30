      module opacity

      use const_lib
      use crlibm_lib
      use chem_def
      use chem_lib
      use eos_def
      use eos_lib

      implicit none

      integer :: handle
      real(dp) :: X, Z, Y, abar, zbar, z2bar, ye
      integer, parameter :: species = 7
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, &
            mg24=7
      integer, pointer, dimension(:) :: net_iso, chem_id
      real(dp) :: xa(species)
      character (len=256) :: my_mesa_dir

!      private
!      public :: init_opacity, shutdown_opacity, get_hbar

      contains

      subroutine init_const()
         integer :: ierr

         ierr = 0
         my_mesa_dir = ''
         call const_init(my_mesa_dir, ierr)
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            write(*,*) my_mesa_dir
            stop 1
         end if
      end subroutine init_const


      subroutine init_chem()
         integer :: ierr

         call crlibm_init
         
         ierr = 0
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in chem_init'
            stop 1
         end if
      end subroutine init_chem


      subroutine init_eos()
         implicit none
         character (len=256) :: eos_file_prefix
         integer :: ierr
         integer, pointer, dimension(:) :: net_iso, chem_id
         logical, parameter :: use_cache = .true.
         integer, parameter :: species = 7
         real(kind=8) :: xa(species)
         real(kind=8), parameter :: X = 0.70
         real(kind=8), parameter :: Z = 0.02
         real(kind=8), parameter :: Zfrac_C = 0.173312d0
         real(kind=8), parameter :: Zfrac_N = 0.053177d0
         real(kind=8), parameter :: Zfrac_O = 0.482398d0
         real(kind=8), parameter :: Zfrac_Ne = 0.098675d0
         real(kind=8) :: frac, dabar_dx(species), dzbar_dx(species), &
               sumx, xh, xhe, xz, mass_correction, dmc_dx(species)

         eos_file_prefix = 'mesa'
         call eos_init(eos_file_prefix, '', '', '', use_cache, ierr)
         if (ierr /= 0) then
            write(*,*) 'eos_init failed'
            stop 1
         end if
         
         handle = alloc_eos_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed trying to allocate eos handle'
            stop 1
         end if
         
         allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
         if (ierr /= 0) stop 'allocate failed'
         net_iso(:) = 0
         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         Y = 1 - (X + Z)               
         xa(h1) = X
         xa(he4) = Y
         xa(c12) = Z * Zfrac_C
         xa(n14) = Z * Zfrac_N
         xa(o16) = Z * Zfrac_O
         xa(ne20) = Z * Zfrac_Ne
         xa(species) = 1 - sum(xa(1:species-1))
         call composition_info( &
               species, chem_id, xa, xh, xhe, xz, abar, zbar, z2bar, &
               ye, mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
      end subroutine init_eos


      subroutine init_opacity()
         implicit none
         call init_const
         call init_chem
         call init_eos
      end subroutine init_opacity


      subroutine shutdown_eos()
         call free_eos_handle(handle)
         call eos_shutdown
         deallocate(net_iso, chem_id)
      end subroutine shutdown_eos


      subroutine shutdown_opacity()
         call shutdown_eos
      end subroutine shutdown_opacity
         


      subroutine get_hbar(x)
         implicit none
         real(kind=8), intent(out) :: x
         
         x = hbar
      end subroutine get_hbar


      end module opacity

