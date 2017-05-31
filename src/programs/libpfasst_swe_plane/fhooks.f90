module hooks_module
  use pf_mod_dtype
  use encap_module
  use feval_module
  implicit none

  interface
     
     ! prototypes of the C functions

     subroutine cecho_error(sd,step) bind(c, name="cecho_error")
       use iso_c_binding
       type(c_ptr),intent(in),value :: sd
       integer,intent(in)           :: step
     end subroutine cecho_error

     subroutine cecho_residual(sd,step,iter) bind(c, name="cecho_residual")
       use iso_c_binding
       type(c_ptr),intent(in),value :: sd
       integer,intent(in)           :: step,iter
     end subroutine cecho_residual
  end interface

contains

  ! error function

  subroutine fecho_error(pf, level, state) 
    use iso_c_binding    
    type(pf_pfasst_t),  intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    type(pf_state_t),   intent(in)    :: state

    class(pf_encap_t),  allocatable   :: Y_exact

    ! allocate memory for vector Y_exact
    call level%ulevel%factory%create_single(Y_exact, & 
                                            level%level, & 
                                            SDC_KIND_SOL_FEVAL, & 
                                            level%nvars, &
                                            level%shape)
    
    ! compute the reference solution (currently with explicit RK)
    call freference(level%ulevel%sweeper, &
                    state%t0+state%dt, & 
                    Y_exact)

    ! compute the norm of the error
    call Y_exact%axpy(-1.0_pfdp, & 
                      level%qend)
    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
         state%step+1, state%iter, Y_exact%norm()

    ! release memory
    call level%ulevel%factory%destroy_single(Y_exact, &
                                             level%level, & 
                                             SDC_KIND_SOL_FEVAL, &
                                             level%nvars, &
                                             level%shape)


  end subroutine fecho_error

  
  ! function to output the residual 

  subroutine fecho_residual(pf, level, state)
    use iso_c_binding 
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state

    select type(R => level%R(level%nnodes-1))
    type is (sweet_data_encap_t)
       
       print '("resid: step: ",i3.3," iter: ",i4.3," level: ",i2.2," resid: ",es14.7)', &
            state%step+1, state%iter, level%level, R%norm()

      !call cecho_residual(R%c_sweet_data_ptr, & 
      !                    state%step+1, & 
      !                    state%iter)

   class default
      stop "TYPE ERROR"
   end select

  end subroutine fecho_residual


end module hooks_module

