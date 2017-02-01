module feval_module
  use iso_c_binding 
  use pf_mod_dtype
  use encap_module
  implicit none

  interface
     subroutine cinitial(ldY, ctx) bind(c, name="cinitial")
       use iso_c_binding
       type(c_ptr),value :: ldY, ctx
     end subroutine cinitial

     subroutine cfinal(ldY, ctx) bind(c, name="cfinal")
       use iso_c_binding
       type(c_ptr),value :: ldY, ctx
     end subroutine cfinal

     subroutine ceval_f1(ldY, t, ldF1, ctx) bind(c, name="ceval_f1")
       use iso_c_binding
       type(c_ptr),value :: ldY, ldF1, ctx
       real(c_double),value :: t
     end subroutine ceval_f1

     subroutine ceval_f2(ldY, t, ldF2, ctx) bind(c, name="ceval_f2")
       use iso_c_binding
       type(c_ptr),value :: ldY, ldF2, ctx
       real(c_double),value :: t
     end subroutine ceval_f2

     subroutine ccomp_f2(ldY, t, dt,ldrhs, ldF2, ctx) bind(c, name="ccomp_f2")
       use iso_c_binding
       type(c_ptr),value :: ldY, ldrhs, ldF2, ctx
       real(c_double),value :: t, dt
     end subroutine ccomp_f2

     subroutine ceval_f3(ldY, t, ldF3, ctx) bind(c, name="ceval_f3")
       use iso_c_binding
       type(c_ptr),value :: ldY, ldF3, ctx
       real(c_double),value :: t
     end subroutine ceval_f3

     subroutine ccomp_f3(ldY, t, dt,ldrhs, ldF3, ctx) bind(c, name="ccomp_f3")
       use iso_c_binding
       type(c_ptr),value :: ldY, ldrhs, ldF3, ctx
       real(c_double),value :: t, dt
     end subroutine ccomp_f3
  end interface

contains

  subroutine finitial(yptr)
    type(c_ptr), intent(in), value :: yptr
    type(level_data), pointer :: y
    call c_f_pointer(yptr,  y)
    call cinitial(y%c_level_data_ptr,y%userCtx_ptr)
  end subroutine finitial

  subroutine ffinal(yptr)
    type(c_ptr), intent(in), value :: yptr
    type(level_data), pointer :: y
    call c_f_pointer(yptr,  y)
    call cfinal(y%c_level_data_ptr,y%userCtx_ptr)
  end subroutine ffinal

  subroutine feval_f1(yptr, t, level, ctx, f1ptr)
    type(c_ptr), intent(in), value :: yptr, f1ptr, ctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level
    type(level_data), pointer :: y, f1

    call c_f_pointer(yptr,  y)
    call c_f_pointer(f1ptr, f1)
    call ceval_f1(y%c_level_data_ptr,t,f1%c_level_data_ptr,y%userCtx_ptr)
  end subroutine feval_f1

  subroutine feval_f2(yptr, t, level, ctx, f2ptr)
    type(c_ptr), intent(in), value :: yptr, f2ptr, ctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level
    type(level_data), pointer :: y, f2
    call c_f_pointer(yptr,  y)
    call c_f_pointer(f2ptr, f2)
    call ceval_f2(y%c_level_data_ptr,t,f2%c_level_data_ptr,y%userCtx_ptr)
  end subroutine feval_f2

  subroutine fcomp_f2(yptr, t, dt, rhsptr, level, ctx, f2ptr)
    type(c_ptr), intent(in), value :: yptr, rhsptr, f2ptr, ctx
    real(pfdp),  intent(in)        :: t, dt
    integer,     intent(in)        :: level
    type(level_data), pointer :: y, rhs, f2
    call c_f_pointer(yptr,      y)
    call c_f_pointer(rhsptr,  rhs)
    call c_f_pointer(f2ptr,    f2)
    call ccomp_f2(y%c_level_data_ptr,t,dt,rhs%c_level_data_ptr,f2%c_level_data_ptr,y%userCtx_ptr)
  end subroutine fcomp_f2

  subroutine feval_f3(yptr, t, level, ctx, f3ptr)
    type(c_ptr), intent(in), value :: yptr, f3ptr, ctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level
    type(level_data), pointer :: y, f3
    call c_f_pointer(yptr,  y)
    call c_f_pointer(f3ptr, f3)
    call ceval_f3(y%c_level_data_ptr,t,f3%c_level_data_ptr,y%userCtx_ptr)
  end subroutine feval_f3

  subroutine fcomp_f3(yptr, t, dt, rhsptr, level, ctx, f3ptr)
    type(c_ptr), intent(in), value :: yptr, rhsptr, f3ptr, ctx
    real(pfdp),  intent(in)        :: t, dt
    integer,     intent(in)        :: level
    type(level_data), pointer :: y, rhs, f3
    call c_f_pointer(yptr,      y)
    call c_f_pointer(rhsptr,  rhs)
    call c_f_pointer(f3ptr,    f3)
    call ccomp_f3(y%c_level_data_ptr,t,dt,rhs%c_level_data_ptr,f3%c_level_data_ptr,y%userCtx_ptr)
  end subroutine fcomp_f3

end module feval_module
