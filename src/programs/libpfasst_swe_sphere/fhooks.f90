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

     subroutine cecho_output_solution(i_ctx, i_Y, i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter) & 
          bind(c, name="cecho_output_solution")
       use iso_c_binding
       type(c_ptr),           value :: i_ctx, i_Y
       integer,               value :: i_current_proc, i_current_step, i_current_iter, i_nnodes, i_niter
     end subroutine cecho_output_solution

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

    select type(R => level%R(level%nnodes-1))
    type is (sweet_data_encap_t)
       
       print '("resid: step: ",i7.5," iter: ",i5.3," level: ",i2.2," resid: ",es14.7)', &
            state%step+1, state%iter, level%level, R%norm()

   class default
      stop "TYPE ERROR"
   end select

  end subroutine fecho_residual


  ! function to output the interpolation errors

  subroutine fecho_interpolation_errors(pf, level, state)
    use iso_c_binding
    use pf_mod_utils
    use pf_mod_restrict
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t), intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state

    class(pf_encap_t), allocatable   :: Q_coarse(:)                  !  Coarse in time and space
    class(pf_encap_t), allocatable   :: Q_fine(:)                    !  Fine in time and space
    class(pf_encap_t), allocatable   :: Q_fine_space_coarse_time(:)  !  Coarse in time but fine in space
    class(pf_encap_t), allocatable   :: Q_coarse_space_fine_time(:)  !  Fine in time but coarse in space

    real(pfdp),        allocatable   :: c_times(:)
    real(pfdp),        allocatable   :: f_times(:)

    integer                          :: m
    logical                          :: interpolate_in_space_first

    interpolate_in_space_first = .true.
    
    if (level%level == pf%nlevels .and. pf%nlevels > 1) then
       
       ! allocate memory for the time nodes
       allocate(c_times(pf%levels(pf%nlevels-1)%nnodes))
       allocate(f_times(pf%levels(pf%nlevels)%nnodes))    
       c_times = state%t0 + state%dt * pf%levels(pf%nlevels-1)%nodes
       f_times = state%t0 + state%dt * pf%levels(pf%nlevels)%nodes
       
       ! allocate memory for the data values
       call pf%levels(pf%nlevels)%ulevel%factory%create_array(Q_fine,                       &
                                                              pf%levels(pf%nlevels)%nnodes, &
                                                              pf%levels(pf%nlevels)%level,  &
                                                              SDC_KIND_CORRECTION,          & 
                                                              pf%levels(pf%nlevels)%nvars,  &
                                                              pf%levels(pf%nlevels)%shape)
       call pf%levels(pf%nlevels-1)%ulevel%factory%create_array(Q_coarse,                       &
                                                                pf%levels(pf%nlevels-1)%nnodes, &
                                                                pf%levels(pf%nlevels-1)%level,  & 
                                                                SDC_KIND_CORRECTION,            &
                                                                pf%levels(pf%nlevels-1)%nvars,  &
                                                                pf%levels(pf%nlevels-1)%shape)
       if (interpolate_in_space_first .eqv. .true.) then
          call pf%levels(pf%nlevels)%ulevel%factory%create_array(Q_fine_space_coarse_time,       &
                                                                 pf%levels(pf%nlevels-1)%nnodes, &
                                                                 pf%levels(pf%nlevels)%level,    &
                                                                 SDC_KIND_CORRECTION,            &
                                                                 pf%levels(pf%nlevels)%nvars,    &
                                                                 pf%levels(pf%nlevels)%shape)
       else
          call pf%levels(pf%nlevels)%ulevel%factory%create_array(Q_coarse_space_fine_time,      &
                                                                 pf%levels(pf%nlevels)%nnodes,  &
                                                                 pf%levels(pf%nlevels-1)%level, &
                                                                 SDC_KIND_CORRECTION,           &
                                                                 pf%levels(pf%nlevels-1)%nvars, &
                                                                 pf%levels(pf%nlevels-1)%shape)
       end if

    
       ! copy the values at all SDC nodes into temporary vectors
       do m = 1, pf%levels(pf%nlevels-1)%nnodes
          call Q_coarse(m)%setval(0.0_pfdp)
       end do
       do m = 1, pf%levels(pf%nlevels)%nnodes
          call Q_fine(m)%copy(pf%levels(pf%nlevels)%Q(m))
       end do
       
       ! restrict in space and then in time
       call restrict_sdc(pf%levels(pf%nlevels),   &
                         pf%levels(pf%nlevels-1), & 
                         Q_fine,                  & 
                         Q_coarse,                &
                         .false.,                 &
                         f_times)
    
       
       do m = 1, pf%levels(pf%nlevels)%nnodes
          call Q_fine(m)%setval(0.0_pfdp)
       end do

       if (interpolate_in_space_first .eqv. .true.) then

          ! interpolate in space
          do m = 1, pf%levels(pf%nlevels-1)%nnodes
             call Q_fine_space_coarse_time(m)%setval(0.0_pfdp)
             call pf%levels(pf%nlevels)%ulevel%interpolate(pf%levels(pf%nlevels),       &
                                                           pf%levels(pf%nlevels-1),     &
                                                           Q_fine_space_coarse_time(m), &
                                                           Q_coarse(m),                 &
                                                           c_times(m))
          end do

          ! interpolate in time
          call pf_apply_mat(Q_fine,                     &
                            1.0_pfdp,                   &
                            pf%levels(pf%nlevels)%tmat, &
                            Q_fine_space_coarse_time,   &
                            .false.)

       else

          ! interpolate in time
          do m = 1, pf%levels(pf%nlevels)%nnodes
             call Q_coarse_space_fine_time(m)%setval(0.0_pfdp)
          end do
          call pf_apply_mat(Q_coarse_space_fine_time,     &
                            1.0_pfdp,                     &
                            pf%levels(pf%nlevels)%tmat,   &
                            Q_coarse,                     &
                            .false.)

          ! interpolate in space
          do m = 1, pf%levels(pf%nlevels)%nnodes
             call Q_fine(m)%setval(0.0_pfdp)
             call pf%levels(pf%nlevels)%ulevel%interpolate(pf%levels(pf%nlevels),       &
                                                           pf%levels(pf%nlevels-1),     &
                                                           Q_fine(m),                   &
                                                           Q_coarse_space_fine_time(m), &
                                                           f_times(m))
          end do


       end if

       ! compute the interpolation error
       do m = 1, pf%levels(pf%nlevels)%nnodes
          call Q_fine(m)%axpy(-1.0_pfdp, pf%levels(pf%nlevels)%Q(m))
          print *, 'interpolation_error: ', Q_fine(m)%norm()
       end do
         
       
       call pf%levels(pf%nlevels-1)%ulevel%factory%destroy_array(Q_coarse,                       &
                                                                 pf%levels(pf%nlevels-1)%nnodes, &
                                                                 pf%levels(pf%nlevels-1)%level,  &
                                                                 SDC_KIND_CORRECTION,            &
                                                                 pf%levels(pf%nlevels-1)%nvars,  &
                                                                 pf%levels(pf%nlevels-1)%shape)
       call pf%levels(pf%nlevels)%ulevel%factory%destroy_array(Q_fine,                       &
                                                               pf%levels(pf%nlevels)%nnodes, &
                                                               pf%levels(pf%nlevels)%level,  &
                                                               SDC_KIND_CORRECTION,          &
                                                               pf%levels(pf%nlevels)%nvars,  &
                                                               pf%levels(pf%nlevels)%shape)

       if (interpolate_in_space_first .eqv. .true.) then
          call pf%levels(pf%nlevels)%ulevel%factory%destroy_array(Q_fine_space_coarse_time,       &
                                                                  pf%levels(pf%nlevels-1)%nnodes, &
                                                                  pf%levels(pf%nlevels)%level,    &
                                                                  SDC_KIND_CORRECTION,            &
                                                                  pf%levels(pf%nlevels)%nvars,    &
                                                                  pf%levels(pf%nlevels)%shape)
       else
          call pf%levels(pf%nlevels)%ulevel%factory%destroy_array(Q_coarse_space_fine_time,       &
                                                                  pf%levels(pf%nlevels)%nnodes,   &
                                                                  pf%levels(pf%nlevels-1)%level,  &
                                                                  SDC_KIND_CORRECTION,            &
                                                                  pf%levels(pf%nlevels-1)%nvars,  &
                                                                  pf%levels(pf%nlevels-1)%shape)
       end if
              
    end if

  end subroutine fecho_interpolation_errors

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

       if (modulo(state%step, 100) == 0 .and. state%iter == pf%niters) then

          call cecho_output_solution(sweet_sweeper_ptr%ctx,  &
                                     x_ptr%c_sweet_data_ptr, &
                                     state%proc,             &
                                     state%step,             &
                                     state%iter,             &
                                     level%nnodes,           &
                                     pf%niters)
          
          print *, 'step = ', state%step, ' iter = ', state%iter, ' processor = ', state%proc

       end if

!    end if

  end subroutine fecho_output_solution

end module hooks_module

