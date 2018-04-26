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

     subroutine cecho_residual(i_ctx, i_norm, i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter) & 
          bind(c, name="cecho_residual")
       use iso_c_binding
       type(c_ptr),     value :: i_ctx
       integer,         value :: i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter
       real(c_double),  value :: i_norm
     end subroutine cecho_residual

     subroutine cecho_output_jump(i_ctx, i_Y, i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter) & 
          bind(c, name="cecho_output_jump")
       use iso_c_binding
       type(c_ptr),     value :: i_ctx, i_Y
       integer,         value :: i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter
     end subroutine cecho_output_jump

     subroutine cecho_output_solution(i_ctx, i_Y, i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter) & 
          bind(c, name="cecho_output_solution")
       use iso_c_binding
       type(c_ptr),     value :: i_ctx, i_Y
       integer,         value :: i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter
     end subroutine cecho_output_solution

     subroutine cecho_output_invariants(i_ctx, i_Y, i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter) & 
          bind(c, name="cecho_output_invariants")
       use iso_c_binding
       type(c_ptr),     value :: i_ctx, i_Y
       integer,         value :: i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter
     end subroutine cecho_output_invariants

  end interface

contains

  ! error function

  subroutine fecho_error(pf, level, state) 
    use iso_c_binding    
    type(pf_pfasst_t),  intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: level
    type(pf_state_t),   intent(in)    :: state

    class(pf_encap_t),  allocatable   :: Y_reference

    ! not implemented yet

  end subroutine fecho_error

  
  ! function to output the residual 

  subroutine fecho_residual(pf, level, state)
    use iso_c_binding 
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: x_ptr

    select type(R => level%R(level%nnodes-1))
    type is (sweet_data_encap_t)
       
       print '("resid: step: ",i7.5," iter: ",i5.3," level: ",i2.2," resid: ",es14.7)', &
            state%step+1, state%iter, level%index, R%norm()

       sweet_sweeper_ptr => as_sweet_sweeper(level%ulevel%sweeper)
       x_ptr             => as_sweet_data_encap(level%Q(sweet_sweeper_ptr%nnodes))

       if (level%index == pf%nlevels) then

          call cecho_residual(sweet_sweeper_ptr%ctx,  &
                              R%norm(),               &
                              state%proc-1,           &
                              state%step,             &
                              state%iter,             &
                              level%nnodes,           &
                              pf%niters)

       end if

   class default
      stop "TYPE ERROR"
   end select

  end subroutine fecho_residual

  ! function to output the jump in the initial condition

  subroutine fecho_output_jump(pf, level, state)
    use iso_c_binding
    use pf_mod_utils
    use pf_mod_restrict
    use mpi
    type(pf_pfasst_t),         intent(inout) :: pf
    class(pf_level_t),         intent(inout) :: level
    type(pf_state_t),          intent(in)    :: state
    
    class(pf_encap_t),         allocatable   :: del
    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: x_ptr
    integer                                  :: ierr, num_procs
    
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

    sweet_sweeper_ptr => as_sweet_sweeper(level%ulevel%sweeper)

    call level%ulevel%factory%create_single(del,                   &
                                            level%index,           & 
                                            SDC_KIND_SOL_NO_FEVAL, &
                                            level%nvars,           &
                                            level%shape)
    call del%copy(level%q0)
    call del%axpy(-1.0_pfdp, level%Q(1))                                                            
    
    x_ptr  => as_sweet_data_encap(del)

    call cecho_output_jump(sweet_sweeper_ptr%ctx,  &
                           x_ptr%c_sweet_data_ptr, &
                           state%proc,             &
                           state%step,             &
                           state%iter,             &
                           level%nnodes,           &
                           pf%niters)

  end subroutine fecho_output_jump
  

  ! function to output the solution

  subroutine fecho_output_solution(pf, level, state)
    use iso_c_binding
    use pf_mod_utils
    use pf_mod_restrict
    use mpi
    type(pf_pfasst_t),         intent(inout) :: pf
    class(pf_level_t),         intent(inout) :: level
    type(pf_state_t),          intent(in)    :: state

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: x_ptr
    integer                                  :: ierr, num_procs

    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
!    if (state%proc == num_procs) then

       sweet_sweeper_ptr => as_sweet_sweeper(level%ulevel%sweeper)
       x_ptr             => as_sweet_data_encap(level%Q(sweet_sweeper_ptr%nnodes))
       
       ! if (modulo(state%step, 10000) == 0 .and. state%iter == pf%niters) then

       call cecho_output_solution(sweet_sweeper_ptr%ctx,  &
                                      x_ptr%c_sweet_data_ptr, &
                                      state%proc,             &
                                      state%step,             &
                                      state%iter,             &
                                      level%nnodes,           &
                                      pf%niters)
          
           !print *, 'step = ', state%step, ' iter = ', state%iter, ' processor = ', state%proc

       ! end if

!    end if

  end subroutine fecho_output_solution

  ! function to output the solution

  subroutine fecho_output_invariants(pf, level, state)
    use iso_c_binding
    use pf_mod_utils
    use pf_mod_restrict
    use mpi
    type(pf_pfasst_t),         intent(inout) :: pf
    class(pf_level_t),         intent(inout) :: level
    type(pf_state_t),          intent(in)    :: state

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: x_ptr
    integer                                  :: ierr, num_procs

    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
!    if (state%proc == num_procs) then

       sweet_sweeper_ptr => as_sweet_sweeper(level%ulevel%sweeper)
       x_ptr             => as_sweet_data_encap(level%Q(sweet_sweeper_ptr%nnodes))

       if (modulo(state%step, 100) == 0 .and. state%iter == pf%niters) then

          call cecho_output_invariants(sweet_sweeper_ptr%ctx,  &
                                       x_ptr%c_sweet_data_ptr, &
                                       state%proc,             &
                                       state%step,             &
                                       state%iter,             &
                                       level%nnodes,           &
                                       pf%niters)

       end if

!    end if
       
     end subroutine fecho_output_invariants


end module hooks_module

