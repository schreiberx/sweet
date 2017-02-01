module hooks_module
  use pf_mod_dtype
  use encap_module
  implicit none

  interface
     subroutine cecho_error(ld,step) bind(c, name="cecho_error")
       use iso_c_binding
       type(c_ptr),intent(in),value :: ld
       integer,intent(in) :: step
     end subroutine cecho_error

     subroutine cecho_residual(ld,step,iter) bind(c, name="cecho_residual")
       use iso_c_binding
       type(c_ptr),intent(in),value :: ld
       integer,intent(in) :: step,iter
     end subroutine cecho_residual
  end interface

contains

  subroutine fecho_error(pf, level, state, ctx)
    use iso_c_binding
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(c_ptr) :: solptr
    type(level_data), pointer :: sol

    solptr = level%Q(1)
    call c_f_pointer(solptr, sol)
    call cecho_error(sol%c_level_data_ptr,state%step+1)
  end subroutine fecho_error

  subroutine echo_residual(pf, level, state, ctx)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(c_ptr) :: resptr
    type(level_data), pointer :: res

    resptr = level%R(level%nnodes-1)
    call c_f_pointer(resptr, res)
    call cecho_residual(res%c_level_data_ptr,state%step+1,state%iter)
  end subroutine echo_residual

end module hooks_module
