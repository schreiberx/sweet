module encap_module
  use iso_c_binding 
  use pf_mod_dtype
  implicit none

  ! derived factory class for the instantion/destruction of the sweet_data_encap_t object

  type, extends (pf_factory_t) :: sweet_data_factory_t
     type(c_ptr) :: ctx = c_null_ptr ! c pointer to PlaneDataCtx/SphereDataCtx
   contains 
     procedure   :: create_single  => sweet_data_create_single
     procedure   :: create_array   => sweet_data_create_array 
     procedure   :: destroy_single => sweet_data_destroy_single
     procedure   :: destroy_array  => sweet_data_destroy_array
  end type sweet_data_factory_t

  type, extends (pf_encap_t) :: sweet_data_encap_t
     type(c_ptr) :: c_sweet_data_ptr = c_null_ptr ! c pointer to PlaneData/SpereData
     integer     :: data_size ! size of the flat data array
   contains
     procedure   :: setval   => sweet_data_setval
     procedure   :: copy     => sweet_data_copy
     procedure   :: norm     => sweet_data_norm
     procedure   :: pack     => sweet_data_pack
     procedure   :: unpack   => sweet_data_unpack
     procedure   :: axpy     => sweet_data_saxpy
     procedure   :: eprint   => sweet_data_eprint
  end type sweet_data_encap_t
  
  interface 

     ! prototypes of the C functions

     subroutine c_sweet_data_create(i_sd_ctx, i_lv, o_sd, o_s) bind(c, name="c_sweet_data_create")
       use iso_c_binding
       type(c_ptr), value       :: i_sd_ctx
       integer,     value       :: i_lv
       type(c_ptr), intent(out) :: o_sd
       integer,     intent(out) :: o_s
     end subroutine c_sweet_data_create

     subroutine c_sweet_data_destroy(io_sd) bind(c, name="c_sweet_data_destroy")
       use iso_c_binding
       type(c_ptr), value :: io_sd
     end subroutine c_sweet_data_destroy

     subroutine c_sweet_data_setval(io_sd, i_val) bind(c, name="c_sweet_data_setval")
       use iso_c_binding
       type(c_ptr),    value :: io_sd
       real(c_double), value :: i_val
     end subroutine c_sweet_data_setval
       
     subroutine c_sweet_data_copy(i_src,o_dst) bind(c, name="c_sweet_data_copy")
       use iso_c_binding
       type(c_ptr), value :: i_src, o_dst
     end subroutine c_sweet_data_copy
     
     subroutine c_sweet_data_norm(i_sd, o_val) bind(c, name="c_sweet_data_norm")
       use iso_c_binding
       type(c_ptr), value :: i_sd
       real(c_double)     :: o_val
     end subroutine c_sweet_data_norm

     subroutine c_sweet_data_pack(io_sd, o_flat_data_ptr) bind(c, name="c_sweet_data_pack")
       use iso_c_binding
       type(c_ptr), value       :: io_sd
       type(c_ptr), intent(out) :: o_flat_data_ptr
     end subroutine c_sweet_data_pack

     subroutine c_sweet_data_unpack(i_flat_data_ptr, o_sd) bind(c, name="c_sweet_data_unpack")
       use iso_c_binding
       type(c_ptr), intent(in) :: i_flat_data_ptr
       type(c_ptr), value      :: o_sd
     end subroutine c_sweet_data_unpack

     subroutine c_sweet_data_saxpy(i_a, i_x, io_y) bind(c, name="c_sweet_data_saxpy")
       use iso_c_binding
       real(c_double), value :: i_a
       type(c_ptr),    value :: i_x, io_y
     end subroutine c_sweet_data_saxpy       

     subroutine c_sweet_data_eprint(i_sd) bind(c, name="c_sweet_data_eprint")
       use iso_c_binding
       type(c_ptr), value :: i_sd
     end subroutine c_sweet_data_eprint

  end interface

contains


  ! function used to "cast" the pf_factory_t objects 
  ! into sweet_data_factory_t objects

  function as_sweet_data_factory(i_factory) result(o_r)
    class(pf_factory_t), intent(in), target :: i_factory
    class(sweet_data_factory_t), pointer    :: o_r

    select type(i_factory)
       type is (sweet_data_factory_t)
          o_r => i_factory
       class default
          stop "TYPE ERROR"
       end select

     end function as_sweet_data_factory


  ! function used to "cast" the pf_encap_t objects 
  ! into sweet_data_encap_t objects
     
  function as_sweet_data_encap(i_encap) result(o_r)
    class(pf_encap_t), intent(in), target :: i_encap
    class(sweet_data_encap_t), pointer    :: o_r

    select type(i_encap)
    type is (sweet_data_encap_t)
       o_r => i_encap
    class default
       stop "TYPE ERROR"
    end select

  end function as_sweet_data_encap
  

  ! constructors/destructors of the sweet_data objects 
  ! PlaneData or SphereData

  subroutine sweet_data_create_single(this, x, level, kind, nvars, shape)
    class(sweet_data_factory_t), intent(inout)              :: this
    class(pf_encap_t),           intent(inout), allocatable :: x
    integer,                     intent(in   )              :: level, kind, nvars, shape(:)

    allocate(sweet_data_encap_t::x)

    select type(x)
    type is (sweet_data_encap_t)
       call c_sweet_data_create(this%ctx, &
                                level-1, &
                                x%c_sweet_data_ptr, &
                                x%data_size) ! conversion to C++ indexing
    class default
       stop "TYPE ERROR"
    end select

  end subroutine sweet_data_create_single


  subroutine sweet_data_create_array(this, x, n, level, kind, nvars, shape)
    class(sweet_data_factory_t), intent(inout)              :: this
    class(pf_encap_t),           intent(inout), allocatable :: x(:)
    integer,                     intent(in   )              :: n, level, kind, nvars, shape(:)

    integer                                                 :: i

    allocate(sweet_data_encap_t::x(n))    

    do i = 1, n
       select type(x(i))
       type is (sweet_data_encap_t)
          call c_sweet_data_create(this%ctx, &
                                   level-1, & 
                                   x(i)%c_sweet_data_ptr, &
                                   x(i)%data_size) ! conversion to C++ indexing    
       class default
          stop "TYPE ERROR"
       end select          
    end do

  end subroutine sweet_data_create_array

  
  subroutine sweet_data_destroy_single(this, x, level, kind, nvars, shape)
    class(sweet_data_factory_t), intent(inout)              :: this
    class(pf_encap_t),           intent(inout), allocatable :: x
    integer,                     intent(in   )              :: level, kind, nvars, shape(:)

    class(sweet_data_encap_t), pointer                      :: x_ptr
    
    x_ptr => as_sweet_data_encap(x)

    call c_sweet_data_destroy(x_ptr%c_sweet_data_ptr)

    deallocate(x)

  end subroutine sweet_data_destroy_single


  subroutine sweet_data_destroy_array(this, x, n, level, kind, nvars, shape)
    class(sweet_data_factory_t), intent(inout)              :: this
    class(pf_encap_t),           intent(inout), allocatable :: x(:)
    integer,                     intent(in   )              :: n, level, kind, nvars, shape(:)

    class(sweet_data_encap_t), pointer                      :: x_ptr
    integer                                                 :: i

    select type(x)
    class is (sweet_data_encap_t)
       do i = 1, n
          call c_sweet_data_destroy(x(i)%c_sweet_data_ptr)
       end do
    end select 
    deallocate(x)

  end subroutine sweet_data_destroy_array


  ! sweet_data (PlaneData/SphereData) operators

  subroutine sweet_data_setval(this, val, flags)
    class(sweet_data_encap_t), intent(inout)           :: this
    real(c_double),            intent(in   )           :: val
    integer,                   intent(in   ), optional :: flags ! not used here

    call c_sweet_data_setval(this%c_sweet_data_ptr, &
                             val)    

  end subroutine sweet_data_setval

  
  subroutine sweet_data_copy(this, src, flags)
    class(sweet_data_encap_t), intent(inout)           :: this
    class(pf_encap_t),         intent(in   )           :: src
    integer,                   intent(in   ), optional :: flags ! not used here

    class(sweet_data_encap_t), pointer                 :: src_sd_ptr

    src_sd_ptr => as_sweet_data_encap(src)

    call c_sweet_data_copy(src_sd_ptr%c_sweet_data_ptr, &
                           this%c_sweet_data_ptr)

  end subroutine sweet_data_copy
    
  
  function sweet_data_norm(this) result (norm)
    class(sweet_data_encap_t), intent(in   ) :: this
    real(c_double)                           :: norm

    call c_sweet_data_norm(this%c_sweet_data_ptr, & 
                           norm)

  end function sweet_data_norm
  

  subroutine sweet_data_saxpy(this, a, x, flags)
    class(sweet_data_encap_t), intent(inout)           :: this
    class(pf_encap_t),         intent(in   )           :: x
    real(pfdp),                intent(in   )           :: a
    integer,                   intent(in   ), optional :: flags ! not used here
  
    class(sweet_data_encap_t), pointer                 :: x_sd_ptr

    x_sd_ptr => as_sweet_data_encap(x)
    
    ! this = i_a * i_x + this
    call c_sweet_data_saxpy(a, & 
                            x_sd_ptr%c_sweet_data_ptr, &
                            this%c_sweet_data_ptr)

  end subroutine sweet_data_saxpy

  
  subroutine sweet_data_pack(this, z)
    class(sweet_data_encap_t), intent(in   ) :: this
    real(pfdp),                intent(  out) :: z(:)

    real(pfdp),                pointer       :: z_ptr(:)
    type(c_ptr)                              :: z_c_ptr

    ! get the data from the sweet object (PlaneData or SphereData)
    call c_sweet_data_pack(this%c_sweet_data_ptr, &
                           z_c_ptr)

    ! convert the C pointer into a Fortran pointer
    call c_f_pointer(z_c_ptr, &
                     z_ptr, &
                     [this%data_size])

    ! copy the data into z
    z = z_ptr

  end subroutine sweet_data_pack


  subroutine sweet_data_unpack(this, z)
    class(sweet_data_encap_t), intent(inout)          :: this
    real(pfdp),                intent(in   ), target  :: z(:)
    real(pfdp),                               target  :: z2(size(z))

    type(c_ptr)                                       :: z_c_ptr

    ! the z2 array is needed because Fortran cannot pass assumed shape arrays to C++
    ! this can probably be optimized
    z2 = z

    ! get a pointer to the first slot in the array
    z_c_ptr = c_loc(z2(1))

    ! copy the data array into the sweet object (PlaneData or SphereData)
    call c_sweet_data_unpack(z_c_ptr, &
                             this%c_sweet_data_ptr)

  end subroutine sweet_data_unpack


  subroutine sweet_data_eprint(this)
    class(sweet_data_encap_t), intent(inout) :: this

    call c_sweet_data_eprint(this%c_sweet_data_ptr)
  end subroutine sweet_data_eprint

end module encap_module

