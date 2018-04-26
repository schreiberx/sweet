module feval_module
  use iso_c_binding 
  use pf_mod_dtype
  use encap_module
  use pf_mod_utils
  use pf_mod_misdcQ
  use pf_mod_rkstepper
  implicit none
  
  ! Define the derived sweeper type
  type, extends(pf_misdcQ_t) :: sweet_sweeper_t
     type(c_ptr)    :: ctx = c_null_ptr ! c pointer to PlaneDataCtx/SphereDataCtx
     integer        :: nnodes           ! number of nodes
     integer        :: sweep_niter      ! number of the current sweep
     integer        :: sweep_niter_max  ! max number of sweeps
     integer        :: level            ! level on which the sweeper is acting 
     real(c_double) :: dt               ! full timestep size
   contains 
     procedure :: f_eval                => sweet_f_eval_sweeper 
     procedure :: f_comp                => sweet_f_comp_sweeper
     procedure :: destroy               => sweet_sweeper_destroy
  end type sweet_sweeper_t

  ! Define the derived sweeper type
  type, extends(pf_ark_t) :: sweet_stepper_t
     type(c_ptr)    :: ctx = c_null_ptr ! c pointer to PlaneDataCtx/SphereDataCtx
     integer        :: nnodes           ! number of nodes
     integer        :: sweep_niter      ! number of the current sweep
     integer        :: sweep_niter_max  ! max number of sweeps
     integer        :: level            ! level on which the sweeper is acting 
     real(c_double) :: dt               ! full timestep size
   contains 
     procedure :: f_eval                => sweet_f_eval_stepper 
     procedure :: f_comp                => sweet_f_comp_stepper
     procedure :: f_finalize            => sweet_f_finalize_stepper
     procedure :: destroy               => sweet_stepper_destroy
  end type sweet_stepper_t
  
  ! prototypes of the C functions

  interface 
     subroutine cinitial(i_ctx, i_t, i_dt, o_Y) bind(c, name="cinitial")
       use iso_c_binding
       type(c_ptr),    value :: i_ctx, o_Y
       real(c_double), value :: i_t, i_dt 
     end subroutine cinitial

     subroutine cfinal(i_ctx, i_Y, i_nnodes, i_niter, i_rank, i_nprocs) bind(c, name="cfinal")
       use iso_c_binding
       type(c_ptr), value :: i_ctx, i_Y
       integer,     value :: i_nnodes
       integer,     value :: i_niter
       integer,     value :: i_rank, i_nprocs
     end subroutine cfinal

     subroutine ceval_f1(i_Y, i_t, i_ctx, o_F1) bind(c, name="ceval_f1")
       use iso_c_binding
       type(c_ptr),    value :: i_Y, i_ctx, o_F1
       real(c_double), value :: i_t
     end subroutine ceval_f1

     subroutine ceval_f2(i_Y, i_t, i_ctx, o_F2) bind(c, name="ceval_f2")
       use iso_c_binding
       type(c_ptr),    value :: i_Y, i_ctx, o_F2
       real(c_double), value :: i_t
     end subroutine ceval_f2

     subroutine ccomp_f2(io_Y, i_t, i_dt, i_Rhs, i_ctx, o_F2) bind(c, name="ccomp_f2")
       use iso_c_binding
       type(c_ptr),    value :: io_Y, i_Rhs, i_ctx, o_F2
       real(c_double), value :: i_t, i_dt
     end subroutine ccomp_f2

     subroutine ceval_f3(i_Y, i_t, i_level, i_ctx, o_F3) bind(c, name="ceval_f3")
       use iso_c_binding
       type(c_ptr),    value :: i_Y, i_ctx, o_F3 
       integer,        value :: i_level
       real(c_double), value :: i_t
     end subroutine ceval_f3

     subroutine ccomp_f3(i_Y, i_t, i_dt, i_level, i_Rhs, i_ctx, o_F3) bind(c, name="ccomp_f3")
       use iso_c_binding
       type(c_ptr),    value :: i_Y, i_Rhs, i_ctx, o_F3
       real(c_double), value :: i_t, i_dt
       integer,        value :: i_level
     end subroutine ccomp_f3

     subroutine cfinalize(i_Y, i_t, i_dt, i_ctx) bind(c, name="cfinalize")
       use iso_c_binding
       type(c_ptr),    value :: i_Y, i_ctx
       real(c_double), value :: i_t, i_dt
     end subroutine cfinalize

  end interface
  
contains

  ! function to cast the pf_sweeper_t object into a sweet_sweeper_t object 
  function as_sweet_sweeper(sweeper) result(r)
    class(pf_sweeper_t),    intent(in), target  :: sweeper
    class(sweet_sweeper_t),             pointer :: r

    select type(sweeper)
    type is (sweet_sweeper_t)
       r => sweeper
    class default
       stop "TYPE ERROR"
    end select

  end function as_sweet_sweeper

  ! function to cast the pf_rkstepper_t object into a sweet_rkstepper_t object
  function as_sweet_stepper(stepper) result(r)
    class(pf_stepper_t),    intent(in), target  :: stepper
    class(sweet_stepper_t),             pointer :: r

    select type(stepper)
    type is (sweet_stepper_t)
       r => stepper
    class default
       stop "TYPE ERROR"
    end select

  end function as_sweet_stepper

  ! initialization of the sweet data vector using z

  subroutine finitial(sweeper, sd, z, t, dt)
    class(pf_sweeper_t),       intent(inout) :: sweeper
    class(pf_encap_t),         intent(inout) :: sd
    class(pf_encap_t),         intent(inout) :: z
    real(c_double),            intent(in)    :: t, dt

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: sd_ptr

    sweet_sweeper_ptr => as_sweet_sweeper(sweeper)
    sd_ptr            => as_sweet_data_encap(sd)

    call cinitial(sweet_sweeper_ptr%ctx, & 
                  t,                     &
                  dt,                    &
                  sd_ptr%c_sweet_data_ptr)

    call z%copy(sd)

  end subroutine finitial


  subroutine ffinal(sweeper, sd, nnodes, niter)
    use mpi 

    class(pf_sweeper_t),       intent(in)    :: sweeper
    class(pf_encap_t),         intent(inout) :: sd
    integer,                   intent(in)    :: nnodes, niter

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: sd_ptr

    integer                                  :: nprocs, rank, ierr

    sweet_sweeper_ptr => as_sweet_sweeper(sweeper)
    sd_ptr            => as_sweet_data_encap(sd)

    call MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs,  ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)


    call cfinal(sweet_sweeper_ptr%ctx,   & 
                sd_ptr%c_sweet_data_ptr, &
                nnodes,                  &
                niter,                   &
                rank,                    &
                nprocs)

  end subroutine ffinal  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!             SWEEPER             !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! evaluate the right-hand side 

  subroutine sweet_f_eval_sweeper(this, y, t, level_index, f, piece)
    class(sweet_sweeper_t),    intent(inout) :: this
    class(pf_encap_t),         intent(in)    :: y
    real(pfdp),                intent(in)    :: t
    integer,                   intent(in)    :: level_index
    class(pf_encap_t),         intent(inout) :: f
    integer,                   intent(in)    :: piece

    class(sweet_data_encap_t), pointer       :: y_sd_ptr
    class(sweet_data_encap_t), pointer       :: f_sd_ptr
    
    y_sd_ptr  => as_sweet_data_encap(y)
    f_sd_ptr  => as_sweet_data_encap(f)

    select case (piece)

       case (1) ! misdc rhs
          call ceval_f1(y_sd_ptr%c_sweet_data_ptr, &
                        t,                         & 
                        this%ctx,                  &  
                        f_sd_ptr%c_sweet_data_ptr)

       case (2) ! first implicit rhs
          call ceval_f2(y_sd_ptr%c_sweet_data_ptr, & 
                        t,                         & 
                        this%ctx,                  & 
                        f_sd_ptr%c_sweet_data_ptr)

       case (3) ! second implicit rhs
          call ceval_f3(y_sd_ptr%c_sweet_data_ptr, & 
                        t,                         & 
                        this%level-1,              &
                        this%ctx,                  &
                        f_sd_ptr%c_sweet_data_ptr)

       case DEFAULT
          print *, 'Piece argument in f_eval can only be 1, 2'!, or 3'
          call exit(0)
    end select    

  end subroutine sweet_f_eval_sweeper

  
  ! solve for y in the implicit system + update f

  subroutine sweet_f_comp_sweeper(this, y, t, dtq, rhs, level_index, f, piece)
    class(sweet_sweeper_t),   intent(inout) :: this
    class(pf_encap_t),        intent(inout) :: y, f
    real(pfdp),               intent(in)    :: t, dtq
    class(pf_encap_t),        intent(in)    :: rhs
    integer,                  intent(in)    :: level_index
    integer,                  intent(in)    :: piece

    class(sweet_data_encap_t), pointer      :: y_sd_ptr
    class(sweet_data_encap_t), pointer      :: f_sd_ptr
    class(sweet_data_encap_t), pointer      :: rhs_sd_ptr

    y_sd_ptr   => as_sweet_data_encap(y)
    f_sd_ptr   => as_sweet_data_encap(f)    
    rhs_sd_ptr => as_sweet_data_encap(rhs) 

    select case (piece)

         case (2) ! first implicit solve
          call ccomp_f2(y_sd_ptr%c_sweet_data_ptr,   & 
                        t,                           & 
                        dtq,                         & 
                        rhs_sd_ptr%c_sweet_data_ptr, &
                        this%ctx,                    & 
                        f_sd_ptr%c_sweet_data_ptr)

        case (3) ! second implicit solve
           call ccomp_f3(y_sd_ptr%c_sweet_data_ptr,   &
                         t,                           & 
                         dtq,                         & 
                         this%level-1,                &
                         rhs_sd_ptr%c_sweet_data_ptr, & 
                         this%ctx,                    & 
                         f_sd_ptr%c_sweet_data_ptr)

       case DEFAULT
          print *, 'Piece argument in f_comp can only be 2'! or 3'
          call exit(0)
    end select
         
  end subroutine sweet_f_comp_sweeper
  
  ! destructor

  subroutine sweet_sweeper_destroy(this, lev)
    class(sweet_sweeper_t), intent(inout) :: this
    class(pf_level_t),      intent(inout) :: lev

    ! need the following line since the "final" keyword is not supported by some (older) compilers
    ! it forces Fortran to destroy the parent class data structures
    call this%misdcQ_destroy(lev) 

  end subroutine sweet_sweeper_destroy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!             STEPPER             !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! evaluate the right-hand side 

  subroutine sweet_f_eval_stepper(this, y, t, level_index, f, piece)
    class(sweet_stepper_t),    intent(inout) :: this
    class(pf_encap_t),         intent(in)    :: y
    real(pfdp),                intent(in)    :: t
    integer,                   intent(in)    :: level_index
    class(pf_encap_t),         intent(inout) :: f
    integer,                   intent(in)    :: piece

    class(sweet_data_encap_t), pointer       :: y_sd_ptr
    class(sweet_data_encap_t), pointer       :: f_sd_ptr
    
    y_sd_ptr  => as_sweet_data_encap(y)
    f_sd_ptr  => as_sweet_data_encap(f)

    select case (piece)

       case (1) ! misdc rhs
          call ceval_f1(y_sd_ptr%c_sweet_data_ptr, &
                        t,                         & 
                        this%ctx,                  &  
                        f_sd_ptr%c_sweet_data_ptr)

       case (2) ! first implicit rhs
          call ceval_f2(y_sd_ptr%c_sweet_data_ptr, & 
                        t,                         & 
                        this%ctx,                  & 
                        f_sd_ptr%c_sweet_data_ptr)

       case (3) ! second implicit rhs
          call ceval_f3(y_sd_ptr%c_sweet_data_ptr, & 
                        t,                         & 
                        level_index-1,             &
                        this%ctx,                  &
                        f_sd_ptr%c_sweet_data_ptr)

       case DEFAULT
          print *, 'Piece argument in f_eval can only be 1, 2'!, or 3'
          call exit(0)
    end select    

  end subroutine sweet_f_eval_stepper

  
  ! solve for y in the implicit system + update f

  subroutine sweet_f_comp_stepper(this, y, t, dtq, rhs, level_index, f, piece)
    class(sweet_stepper_t),   intent(inout) :: this
    class(pf_encap_t),        intent(inout) :: y, f
    real(pfdp),               intent(in)    :: t, dtq
    class(pf_encap_t),        intent(in)    :: rhs
    integer,                  intent(in)    :: level_index
    integer,                  intent(in)    :: piece

    class(sweet_data_encap_t), pointer      :: y_sd_ptr
    class(sweet_data_encap_t), pointer      :: f_sd_ptr
    class(sweet_data_encap_t), pointer      :: rhs_sd_ptr

    y_sd_ptr   => as_sweet_data_encap(y)
    f_sd_ptr   => as_sweet_data_encap(f)    
    rhs_sd_ptr => as_sweet_data_encap(rhs) 

    select case (piece)

         case (2) ! first implicit solve
          call ccomp_f2(y_sd_ptr%c_sweet_data_ptr,   & 
                        t,                           & 
                        dtq,                          & 
                        rhs_sd_ptr%c_sweet_data_ptr, &
                        this%ctx,                    & 
                        f_sd_ptr%c_sweet_data_ptr)

        case (3) ! second implicit solve
           call ccomp_f3(y_sd_ptr%c_sweet_data_ptr,   &
                         t,                           & 
                         dtq,                          & 
                         level_index-1,               &
                         rhs_sd_ptr%c_sweet_data_ptr, & 
                         this%ctx,                    & 
                         f_sd_ptr%c_sweet_data_ptr)

       case DEFAULT
          print *, 'Piece argument in f_comp can only be 2'! or 3'
          call exit(0)
    end select
         
  end subroutine sweet_f_comp_stepper

  ! apply artificial viscosity
  
    subroutine sweet_f_finalize_stepper(this, y, t, dtq, level_index)
    class(sweet_stepper_t), intent(inout) :: this
    class(pf_encap_t),      intent(inout) :: y
    real(pfdp),             intent(in   ) :: t
    real(pfdp),             intent(in   ) :: dtq
    integer,                intent(in   ) :: level_index

    class(sweet_data_encap_t),    pointer :: y_sd_ptr

    call cfinalize(y_sd_ptr%c_sweet_data_ptr, &
                   t,                         &
                   dtq,                       &
                   this%ctx)                   
    
  end subroutine sweet_f_finalize_stepper
  
  ! destructor

  subroutine sweet_stepper_destroy(this, lev)
    class(sweet_stepper_t), intent(inout) :: this
    class(pf_level_t),      intent(inout) :: lev

    ! need the following line since the "final" keyword is not supported by some (older) compilers
    ! it forces Fortran to destroy the parent class data structures
    call this%ark_destroy(lev) 

  end subroutine sweet_stepper_destroy


end module feval_module

