module encap_module

  use iso_c_binding 
  use pf_mod_dtype
  implicit none

  type :: level_data
     type(c_ptr) :: c_level_data_ptr
     type(c_ptr) :: userCtx_ptr
     integer :: ncomp, level
  end type level_data

  interface
     subroutine c_level_data_create(ctx,ld,lv) bind(c, name="c_level_data_create")
       use iso_c_binding
       type(c_ptr),value :: ctx
       type(c_ptr),intent(inout) :: ld
       integer,value:: lv
     end subroutine c_level_data_create

     subroutine c_level_data_destroy(ld) bind(c, name="c_level_data_destroy")
       use iso_c_binding
       type(c_ptr),value :: ld
     end subroutine c_level_data_destroy

     subroutine c_level_data_setval(ld,val) bind(c, name="c_level_data_setval")
       use iso_c_binding
       type(c_ptr),value:: ld
       real(c_double),value:: val
     end subroutine c_level_data_setval

     subroutine c_level_data_copy(dst,src) bind(c, name="c_level_data_copy")
       use iso_c_binding
       type(c_ptr),value:: dst, src
     end subroutine c_level_data_copy

     subroutine c_level_data_norm(ld,val) bind(c, name="c_level_data_norm")
       use iso_c_binding
       type(c_ptr),value:: ld
       real(c_double) :: val
     end subroutine c_level_data_norm

     subroutine c_level_data_saxpy(y,a,x) bind(c, name="c_level_data_saxpy")
       use iso_c_binding
       type(c_ptr),value:: y,x
       real(c_double),value :: a
     end subroutine c_level_data_saxpy

     subroutine c_level_data_dump(ld, name) bind(c, name="c_level_data_dump")
       use iso_c_binding
       character(c_char) :: name(*)
       type(c_ptr),value :: ld
     end subroutine c_level_data_dump

  end interface

contains

  subroutine c_encap_create(solptr, level, kind, nvars, shape, ctx)
    type(c_ptr),      intent(inout) :: solptr
    integer,          intent(in   ) :: level, nvars, shape(:)
    integer,          intent(in   ) :: kind
    type(c_ptr),      intent(in   ),value :: ctx
    type(level_data), pointer       :: sol
    allocate(sol)
    sol%ncomp = nvars
    sol%level = level - 1
    sol%userCtx_ptr = ctx
    call c_level_data_create(sol%userCtx_ptr,sol%c_level_data_ptr,sol%level)
    solptr = c_loc(sol)
  end subroutine c_encap_create

  subroutine c_encap_destroy(xptr)
    type(c_ptr), intent(in   ), value :: xptr
    type(level_data), pointer :: x
    call c_f_pointer(xptr, x)
    call c_level_data_destroy(x%c_level_data_ptr)
    deallocate(x)
  end subroutine c_encap_destroy

  subroutine c_encap_setval(solptr, val, flags)
    type(c_ptr),      intent(in   ), value    :: solptr
    real(c_double),   intent(in   )           :: val
    integer,          intent(in   ), optional :: flags
    type(level_data), pointer                 :: sol
    call c_f_pointer(solptr, sol)
    call c_level_data_setval(sol%c_level_data_ptr,val)    
  end subroutine c_encap_setval

  subroutine c_encap_copy(dstptr, srcptr, flags)
    type(c_ptr), intent(in   ), value    :: dstptr, srcptr
    integer,     intent(in   ), optional :: flags
    type(level_data), pointer            :: dst, src
    call c_f_pointer(dstptr, dst)
    call c_f_pointer(srcptr, src)
    call c_level_data_copy(dst%c_level_data_ptr,src%c_level_data_ptr)
  end subroutine c_encap_copy

  function c_encap_norm(ptr) result (norm)
    type(c_ptr), intent(in   ), value :: ptr
    real(c_double)                    :: norm
    type(level_data), pointer         :: y
    call c_f_pointer(ptr, y)
    call c_level_data_norm(y%c_level_data_ptr,norm)
  end function c_encap_norm

  subroutine c_encap_saxpy(yptr, a, xptr, flags)
    type(c_ptr), intent(in   ), value    :: yptr, xptr
    real(pfdp), intent(in)               :: a 
    integer,     intent(in   ), optional :: flags
    type(level_data), pointer            :: y, x
    call c_f_pointer(yptr, y)
    call c_f_pointer(xptr, x)
    call c_level_data_saxpy(y%c_level_data_ptr,a,x%c_level_data_ptr)
  end subroutine c_encap_saxpy

  subroutine c_encap_pack(z, ptr)
    type(c_ptr), intent(in   ), value :: ptr
    real(pfdp),  intent(  out)        :: z(:)
    ! do nothing
    print*,'Pack called'
  end subroutine c_encap_pack

  subroutine c_encap_unpack(ptr, z)
    type(c_ptr), intent(in   ), value :: ptr
    real(pfdp),  intent(in   )        :: z(:)
    ! do nothing
    print*,'Unpack called'
  end subroutine c_encap_unpack

  subroutine c_encap_dump(xptr, name, nl)
    type(c_ptr), intent(in   ), value :: xptr
    integer, intent(in) :: nl
    character(c_char) :: name(nl)
    type(level_data), pointer :: x
    character*72 :: tmp
    integer :: i
    do i=1,nl
       tmp(i:i) = name(i)
    enddo
    call c_f_pointer(xptr, x)
    call c_level_data_dump(x%c_level_data_ptr, tmp//c_null_char)
  end subroutine c_encap_dump

  subroutine c_encap_encap_create(encap)
    type(pf_encap_t), intent(out) :: encap
    encap%create  => c_encap_create
    encap%destroy => c_encap_destroy
    encap%setval  => c_encap_setval
    encap%copy    => c_encap_copy
    encap%norm    => c_encap_norm
    encap%pack    => c_encap_pack
    encap%unpack  => c_encap_unpack
    encap%axpy    => c_encap_saxpy
  end subroutine c_encap_encap_create
end module encap_module
