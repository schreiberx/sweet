
module transfer_module
  use pf_mod_dtype
  use encap_module
  use feval_module
  implicit none

  interface
     subroutine c_level_data_restrict(ldG, ldF, t) bind(c, name="c_level_data_restrict")
       use iso_c_binding
       type(c_ptr),value:: ldG, ldF
       double precision,value :: t
     end subroutine c_level_data_restrict

     subroutine c_level_data_interpolate(ldF, ldG, t) bind(c, name="c_level_data_interpolate")
       use iso_c_binding
       type(c_ptr),value:: ldF, ldG
       double precision,value :: t
     end subroutine c_level_data_interpolate
  end interface
  
contains

  subroutine c_encap_interpolate(qFp, qGp, levelF, ctxF, levelG, ctxG, t)
    type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG
    real(pfdp),  intent(in)        :: t
    type(level_data), pointer      :: qF, qG
    call c_f_pointer(qFp, qF)
    call c_f_pointer(qGp, qG)
    call c_level_data_interpolate(qF%c_level_data_ptr, qG%c_level_data_ptr, t)
  end subroutine c_encap_interpolate

  subroutine c_encap_restrict(qFp, qGp, levelF, ctxF, levelG, ctxG, t)
    type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG
    real(pfdp),  intent(in)        :: t
    type(level_data), pointer      :: qF, qG
    call c_f_pointer(qFp, qF)
    call c_f_pointer(qGp, qG)
    call c_level_data_restrict(qG%c_level_data_ptr, qF%c_level_data_ptr, t)
  end subroutine c_encap_restrict
end module transfer_module
