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

     subroutine cecho_residual(i_ctx, i_norm, i_current_proc) &
          bind(c, name="cecho_residual")
       use iso_c_binding
       type(c_ptr),     value :: i_ctx
       integer,         value :: i_current_proc
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

  subroutine fecho_error(pf, level_index)
    use iso_c_binding
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in)              :: level_index

    class(pf_encap_t),  allocatable   :: Y_reference

    ! not implemented yet

  end subroutine fecho_error

  
  ! function to output the residual 

  subroutine fecho_residual(pf, level_index)
    use iso_c_binding 
    use pf_mod_utils
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in)              :: level_index

    real(pfdp)                       :: resid
    integer                          :: step,rank,iter

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: x_ptr

    step=pf%state%step+1
    rank=pf%rank
    iter=pf%state%iter
    resid=pf%levels(level_index)%residual

    print '("[MULE] libpfasst.resid_s",i5.5,"_i",i3.3,"_l",i2.2," = ",es14.7)', &
                step, iter, level_index, resid

    sweet_sweeper_ptr => as_sweet_sweeper(pf%levels(level_index)%ulevel%sweeper)
    x_ptr             => as_sweet_data_encap(pf%levels(level_index)%Q(sweet_sweeper_ptr%nnodes))

   if (level_index == pf%nlevels) then

      call cecho_residual(sweet_sweeper_ptr%ctx,  &
                          resid,               &
                          pf%state%proc-1)

   end if

  end subroutine fecho_residual

  ! function to output the jump in the initial condition

  subroutine fecho_output_jump(pf, level_index)
    use iso_c_binding
    use pf_mod_utils
    use pf_mod_restrict
    use mpi
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in)              :: level_index
    
    class(pf_encap_t),         allocatable   :: del
    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: x_ptr
    integer                                  :: ierr, num_procs

    integer                          ::   proc,step,rank,iter
    
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

    sweet_sweeper_ptr => as_sweet_sweeper(pf%levels(level_index)%ulevel%sweeper)

    proc=pf%state%proc
    step=pf%state%step
    rank=pf%rank
    iter=pf%state%iter

    call pf%levels(level_index)%ulevel%factory%create_single(del,  &
                                            level_index,           &
                                            pf%levels(level_index)%lev_shape)
    call del%copy(pf%levels(level_index)%q0)
    call del%axpy(-1.0_pfdp, pf%levels(level_index)%Q(1))
    
    x_ptr  => as_sweet_data_encap(del)

    call cecho_output_jump(sweet_sweeper_ptr%ctx,  &
                           x_ptr%c_sweet_data_ptr, &
                           proc,                   &
                           step,                   &
                           iter,                   &
                           pf%levels(level_index)%nnodes,           &
                           pf%niters)

  end subroutine fecho_output_jump
  

  ! function to output the solution

  subroutine fecho_output_solution(pf, level_index)
    use iso_c_binding
    use pf_mod_utils
    use pf_mod_restrict
    use mpi
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in)              :: level_index

    integer                          ::   proc,step,rank,iter

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: x_ptr
    integer                                  :: ierr, num_procs

    proc=pf%state%proc
    step=pf%state%step
    rank=pf%rank
    iter=pf%state%iter

    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

       sweet_sweeper_ptr => as_sweet_sweeper(pf%levels(level_index)%ulevel%sweeper)
       x_ptr             => as_sweet_data_encap(pf%levels(level_index)%Q(sweet_sweeper_ptr%nnodes))

       call cecho_output_solution(sweet_sweeper_ptr%ctx,  &
                                      x_ptr%c_sweet_data_ptr, &
                                      proc,             &
                                      step,             &
                                      iter,             &
                                      pf%levels(level_index)%nnodes,           &
                                      pf%niters)

  end subroutine fecho_output_solution

  ! function to output the solution

  subroutine fecho_output_invariants(pf, level_index)
    use iso_c_binding
    use pf_mod_utils
    use pf_mod_restrict
    use mpi
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in)              :: level_index

    integer                          :: proc,step,rank,iter

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: x_ptr
    integer                                  :: ierr, num_procs

    proc=pf%state%proc
    step=pf%state%step
    rank=pf%rank
    iter=pf%state%iter

    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
!    if (state%proc == num_procs) then

       sweet_sweeper_ptr => as_sweet_sweeper(pf%levels(level_index)%ulevel%sweeper)
       x_ptr             => as_sweet_data_encap(pf%levels(level_index)%Q(sweet_sweeper_ptr%nnodes))

       if (modulo(step, 100) == 0 .and. iter == pf%niters) then

          call cecho_output_invariants(sweet_sweeper_ptr%ctx,  &
                                       x_ptr%c_sweet_data_ptr, &
                                       proc,             &
                                       step,             &
                                       iter,             &
                                       pf%levels(level_index)%nnodes,           &
                                       pf%niters)

       end if
       
     end subroutine fecho_output_invariants

end module hooks_module

