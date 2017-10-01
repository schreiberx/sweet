module feval_module
  use iso_c_binding 
  use pf_mod_dtype
  use encap_module
  use pf_mod_utils
  use pf_mod_misdcQ
  implicit none
  
  type, extends(pf_misdcQ_t) :: sweet_sweeper_t
     type(c_ptr)    :: ctx = c_null_ptr ! c pointer to PlaneDataCtx/SphereDataCtx
     integer        :: nnodes           ! number of nodes
     integer        :: sweep_niter      ! number of the current sweep
     integer        :: sweep_niter_max  ! max number of sweeps
     integer        :: level            ! level on which the sweeper is acting 
     real(c_double) :: dt               ! full timestep size
   contains 
     procedure :: f_eval                => sweet_f_eval 
     procedure :: f_comp                => sweet_f_comp
     procedure :: initialize_correction => sweet_initialize_correction
     procedure :: finalize_correction   => sweet_finalize_correction
     procedure :: destroy               => sweet_sweeper_destroy
  end type sweet_sweeper_t
  
  ! prototypes of the C functions

  interface 
     subroutine cinitial(i_ctx, i_t, i_dt, o_Y) bind(c, name="cinitial")
       use iso_c_binding
       type(c_ptr),    value :: i_ctx, o_Y
       real(c_double), value :: i_t, i_dt 
     end subroutine cinitial

     subroutine cfinal(i_ctx, i_Y, i_nnodes, i_niter) bind(c, name="cfinal")
       use iso_c_binding
       type(c_ptr), value :: i_ctx, i_Y
       integer,     value :: i_nnodes
       integer,     value :: i_niter
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

     subroutine capply_viscosity(io_Y, i_t, i_dt, i_nnodes, i_level, i_ctx) bind(c, name="capply_viscosity")
       use iso_c_binding
       type(c_ptr),    value :: io_Y, i_ctx
       real(c_double), value :: i_t, i_dt
       integer,        value :: i_level, i_nnodes
     end subroutine capply_viscosity

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

  end subroutine finitial


  subroutine ffinal(sweeper, sd, nnodes, niter)
    class(pf_sweeper_t),       intent(in)    :: sweeper
    class(pf_encap_t),         intent(inout) :: sd
    integer,                   intent(in)    :: nnodes, niter

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: sd_ptr

    sweet_sweeper_ptr => as_sweet_sweeper(sweeper)
    sd_ptr            => as_sweet_data_encap(sd)

    call cfinal(sweet_sweeper_ptr%ctx,   & 
                sd_ptr%c_sweet_data_ptr, &
                nnodes,                  &
                niter)

  end subroutine ffinal


  ! prepare the current solution at m for the correction                                                                                                                           
  subroutine sweet_initialize_correction(this, y, t, dt, level, m, flag)
    class(sweet_sweeper_t), intent(inout) :: this
    class(pf_encap_t),      intent(inout) :: y
    real(pfdp),             intent(in)    :: t, dt
    integer,                intent(in)    :: level, m
    logical,                intent(in)    :: flag

    ! not implemented

  end subroutine sweet_initialize_correction
  

  ! evaluate the right-hand side 

  subroutine sweet_f_eval(this, y, t, level, m, f, piece)
    class(sweet_sweeper_t),    intent(inout) :: this
    class(pf_encap_t),         intent(in)    :: y
    real(pfdp),                intent(in)    :: t
    integer,                   intent(in)    :: level, m
    class(pf_encap_t),         intent(inout) :: f
    integer,                   intent(in)    :: piece

    class(sweet_data_encap_t), pointer       :: y_sd_ptr
    class(sweet_data_encap_t), pointer       :: f_sd_ptr
    
    y_sd_ptr  => as_sweet_data_encap(y)
    f_sd_ptr  => as_sweet_data_encap(f)

    select case (piece)

       case (1) ! imex rhs
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

  end subroutine sweet_f_eval

  
  ! solve for y in the implicit system + update f

  subroutine sweet_f_comp(this, y, t, dt, rhs, level, f, piece)
    class(sweet_sweeper_t),   intent(inout) :: this
    class(pf_encap_t),        intent(inout) :: y, f
    real(pfdp),               intent(in)    :: t, dt
    class(pf_encap_t),        intent(in)    :: rhs
    integer,                  intent(in)    :: level
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
                        dt,                          & 
                        rhs_sd_ptr%c_sweet_data_ptr, &
                        this%ctx,                    & 
                        f_sd_ptr%c_sweet_data_ptr)

        case (3) ! second implicit solve
           call ccomp_f3(y_sd_ptr%c_sweet_data_ptr,   &
                         t,                           & 
                         dt,                          & 
                         this%level-1,                &
                         rhs_sd_ptr%c_sweet_data_ptr, & 
                         this%ctx,                    & 
                         f_sd_ptr%c_sweet_data_ptr)

       case DEFAULT
          print *, 'Piece argument in f_comp can only be 2'! or 3'
          call exit(0)
    end select
         
  end subroutine sweet_f_comp

 ! apply various treatments to the solution after the correction                                                                                                                  
  subroutine sweet_finalize_correction(this, y, t, dt, level, m, flag)
    class(sweet_sweeper_t), intent(inout) :: this
    class(pf_encap_t),      intent(inout) :: y
    real(pfdp),             intent(in)    :: t, dt
    integer,                intent(in)    :: level, m
    logical,                intent(in)    :: flag 
    
    class(sweet_data_encap_t), pointer    :: y_sd_ptr

    ! if m is the last node to be swept                                                                                                                                         
    if (m == this%nnodes) then

       ! increment the sweep number                                                                                                                                             
       this%sweep_niter = this%sweep_niter + 1

       ! set the delta_chi zero when we start a new timestep                                                                                                                    
       if (this%sweep_niter == this%sweep_niter_max) then
          
          this%sweep_niter = 0 ! reset the counter

          ! do nothing for now

       end if
       
    endif

  end subroutine sweet_finalize_correction
  
  ! destructor

  subroutine sweet_sweeper_destroy(this, lev)
    class(sweet_sweeper_t), intent(inout) :: this
    class(pf_level_t),      intent(inout) :: lev

    ! need the following line since the "final" keyword is not supported by some (older) compilers
    ! it forces Fortran to destroy the parent class data structures
    call this%misdcQ_destroy(lev) 

  end subroutine sweet_sweeper_destroy

end module feval_module

