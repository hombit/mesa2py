! Module class_opacity defined in file opacity.f90

subroutine f90wrap_opacity__get__handle(this, f90wrap_handle)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_handle
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_handle = this_ptr%p%handle
end subroutine f90wrap_opacity__get__handle

subroutine f90wrap_opacity__set__handle(this, f90wrap_handle)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_handle
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%handle = f90wrap_handle
end subroutine f90wrap_opacity__set__handle

subroutine f90wrap_opacity__array__xa(this, nd, dtype, dshape, dloc)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in) :: this(2)
    type(opacity_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%xa)
    dloc = loc(this_ptr%p%xa)
end subroutine f90wrap_opacity__array__xa

subroutine f90wrap_opacity__get__y(this, f90wrap_y)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y = this_ptr%p%y
end subroutine f90wrap_opacity__get__y

subroutine f90wrap_opacity__set__y(this, f90wrap_y)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y = f90wrap_y
end subroutine f90wrap_opacity__set__y

subroutine f90wrap_opacity__get__abar(this, f90wrap_abar)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_abar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_abar = this_ptr%p%abar
end subroutine f90wrap_opacity__get__abar

subroutine f90wrap_opacity__set__abar(this, f90wrap_abar)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_abar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%abar = f90wrap_abar
end subroutine f90wrap_opacity__set__abar

subroutine f90wrap_opacity__get__zbar(this, f90wrap_zbar)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_zbar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zbar = this_ptr%p%zbar
end subroutine f90wrap_opacity__get__zbar

subroutine f90wrap_opacity__set__zbar(this, f90wrap_zbar)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_zbar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%zbar = f90wrap_zbar
end subroutine f90wrap_opacity__set__zbar

subroutine f90wrap_opacity__get__z2bar(this, f90wrap_z2bar)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_z2bar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_z2bar = this_ptr%p%z2bar
end subroutine f90wrap_opacity__get__z2bar

subroutine f90wrap_opacity__set__z2bar(this, f90wrap_z2bar)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_z2bar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%z2bar = f90wrap_z2bar
end subroutine f90wrap_opacity__set__z2bar

subroutine f90wrap_opacity__get__ye(this, f90wrap_ye)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_ye
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ye = this_ptr%p%ye
end subroutine f90wrap_opacity__get__ye

subroutine f90wrap_opacity__set__ye(this, f90wrap_ye)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_ye
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ye = f90wrap_ye
end subroutine f90wrap_opacity__set__ye

subroutine f90wrap_opacity__get__x(this, f90wrap_x)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x = this_ptr%p%x
end subroutine f90wrap_opacity__get__x

subroutine f90wrap_opacity__get__z(this, f90wrap_z)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_z
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_z = this_ptr%p%z
end subroutine f90wrap_opacity__get__z

subroutine f90wrap_opacity__get__zfrac_c(this, f90wrap_zfrac_c)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_zfrac_c
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zfrac_c = this_ptr%p%zfrac_c
end subroutine f90wrap_opacity__get__zfrac_c

subroutine f90wrap_opacity__get__zfrac_n(this, f90wrap_zfrac_n)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_zfrac_n
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zfrac_n = this_ptr%p%zfrac_n
end subroutine f90wrap_opacity__get__zfrac_n

subroutine f90wrap_opacity__get__zfrac_o(this, f90wrap_zfrac_o)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_zfrac_o
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zfrac_o = this_ptr%p%zfrac_o
end subroutine f90wrap_opacity__get__zfrac_o

subroutine f90wrap_opacity__get__zfrac_ne(this, f90wrap_zfrac_ne)
    use class_opacity, only: opacity
    implicit none
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    integer, intent(in)   :: this(2)
    type(opacity_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_zfrac_ne
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zfrac_ne = this_ptr%p%zfrac_ne
end subroutine f90wrap_opacity__get__zfrac_ne

subroutine f90wrap_opacity_initialise(this)
    use class_opacity, only: opacity
    implicit none
    
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    type(opacity_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_opacity_initialise

subroutine f90wrap_opacity_finalise(this)
    use class_opacity, only: opacity
    implicit none
    
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    type(opacity_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_opacity_finalise

subroutine f90wrap_init_const
    use class_opacity, only: init_const
    implicit none
    
    call init_const()
end subroutine f90wrap_init_const

subroutine f90wrap_init_chem
    use class_opacity, only: init_chem
    implicit none
    
    call init_chem()
end subroutine f90wrap_init_chem

subroutine f90wrap_init_eos(op)
    use class_opacity, only: init_eos, opacity
    implicit none
    
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    type(opacity_ptr_type) :: op_ptr
    integer, intent(in), dimension(2) :: op
    op_ptr = transfer(op, op_ptr)
    call init_eos(op=op_ptr%p)
end subroutine f90wrap_init_eos

subroutine f90wrap_init_opacity(ret_init_opacity)
    use class_opacity, only: init_opacity, opacity
    implicit none
    
    type opacity_ptr_type
        type(opacity), pointer :: p => NULL()
    end type opacity_ptr_type
    type(opacity_ptr_type) :: ret_init_opacity_ptr
    integer, intent(out), dimension(2) :: ret_init_opacity
    allocate(ret_init_opacity_ptr%p)
    ret_init_opacity_ptr%p = init_opacity()
    ret_init_opacity = transfer(ret_init_opacity_ptr, ret_init_opacity)
end subroutine f90wrap_init_opacity

subroutine f90wrap_shutdown_eos
    use class_opacity, only: shutdown_eos
    implicit none
    
    call shutdown_eos()
end subroutine f90wrap_shutdown_eos

subroutine f90wrap_shutdown_opacity
    use class_opacity, only: shutdown_opacity
    implicit none
    
    call shutdown_opacity()
end subroutine f90wrap_shutdown_opacity

subroutine f90wrap_get_hbar(x)
    use class_opacity, only: get_hbar
    implicit none
    
    real(4), intent(out) :: x
    call get_hbar(x=x)
end subroutine f90wrap_get_hbar

subroutine f90wrap_class_opacity__get__use_cache(f90wrap_use_cache)
    use class_opacity, only: class_opacity_use_cache => use_cache
    implicit none
    logical, intent(out) :: f90wrap_use_cache
    
    f90wrap_use_cache = class_opacity_use_cache
end subroutine f90wrap_class_opacity__get__use_cache

subroutine f90wrap_class_opacity__get__species(f90wrap_species)
    use class_opacity, only: class_opacity_species => species
    implicit none
    integer, intent(out) :: f90wrap_species
    
    f90wrap_species = class_opacity_species
end subroutine f90wrap_class_opacity__get__species

subroutine f90wrap_class_opacity__get__h1(f90wrap_h1)
    use class_opacity, only: class_opacity_h1 => h1
    implicit none
    integer, intent(out) :: f90wrap_h1
    
    f90wrap_h1 = class_opacity_h1
end subroutine f90wrap_class_opacity__get__h1

subroutine f90wrap_class_opacity__get__he4(f90wrap_he4)
    use class_opacity, only: class_opacity_he4 => he4
    implicit none
    integer, intent(out) :: f90wrap_he4
    
    f90wrap_he4 = class_opacity_he4
end subroutine f90wrap_class_opacity__get__he4

subroutine f90wrap_class_opacity__get__c12(f90wrap_c12)
    use class_opacity, only: class_opacity_c12 => c12
    implicit none
    integer, intent(out) :: f90wrap_c12
    
    f90wrap_c12 = class_opacity_c12
end subroutine f90wrap_class_opacity__get__c12

subroutine f90wrap_class_opacity__get__n14(f90wrap_n14)
    use class_opacity, only: class_opacity_n14 => n14
    implicit none
    integer, intent(out) :: f90wrap_n14
    
    f90wrap_n14 = class_opacity_n14
end subroutine f90wrap_class_opacity__get__n14

subroutine f90wrap_class_opacity__get__o16(f90wrap_o16)
    use class_opacity, only: class_opacity_o16 => o16
    implicit none
    integer, intent(out) :: f90wrap_o16
    
    f90wrap_o16 = class_opacity_o16
end subroutine f90wrap_class_opacity__get__o16

subroutine f90wrap_class_opacity__get__ne20(f90wrap_ne20)
    use class_opacity, only: class_opacity_ne20 => ne20
    implicit none
    integer, intent(out) :: f90wrap_ne20
    
    f90wrap_ne20 = class_opacity_ne20
end subroutine f90wrap_class_opacity__get__ne20

subroutine f90wrap_class_opacity__get__mg24(f90wrap_mg24)
    use class_opacity, only: class_opacity_mg24 => mg24
    implicit none
    integer, intent(out) :: f90wrap_mg24
    
    f90wrap_mg24 = class_opacity_mg24
end subroutine f90wrap_class_opacity__get__mg24

! End of module class_opacity defined in file opacity.f90

