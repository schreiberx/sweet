module feval_module
    use iso_c_binding 
    use pf_mod_dtype
    use encap_module
    use pf_mod_utils
    use pf_mod_imex_sweeper
    use pf_mod_rkstepper
    implicit none
  
    ! Define the derived sweeper type
    type, extends(pf_imex_sweeper_t) :: sweet_sweeper_t
        type(c_ptr)    :: ctx = c_null_ptr ! c pointer to PlaneDataCtx/SphereDataCtx
        integer        :: nnodes           ! number of nodes
        integer        :: sweep_niter      ! number of the current sweep
        integer        :: sweep_niter_max  ! max number of sweeps
        integer        :: level            ! level on which the sweeper is acting 
        real(c_double) :: dt               ! full timestep size
    contains 
        procedure :: f_eval                => sweet_f_eval
        procedure :: f_comp                => sweet_f_comp
        procedure :: initialize            => sweet_sweeper_initialize
        procedure :: destroy               => sweet_sweeper_destroy
        procedure :: compute_dt            => sweet_sweeper_compute_dt
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

        subroutine ceval(i_Y, i_t, i_ctx, o_F) bind(c, name="ceval")
            use iso_c_binding
            type(c_ptr),    value :: i_Y, i_ctx, o_F
            real(c_double), value :: i_t
        end subroutine ceval

        subroutine ccomp(io_Y, i_t, i_dt, i_Rhs, i_ctx, o_F) bind(c, name="ccomp")
            use iso_c_binding
            type(c_ptr),    value :: io_Y, i_Rhs, i_ctx, o_F
            real(c_double), value :: i_t, i_dt
        end subroutine ccomp

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

        sweet_sweeper_ptr => as_sweet_sweeper(sweeper)
        sd_ptr            => as_sweet_data_encap(sd)

        call cfinal(sweet_sweeper_ptr%ctx,   & 
                    sd_ptr%c_sweet_data_ptr, &
                    nnodes,                  &
                    niter)

    end subroutine ffinal  

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!             SWEEPER             !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! evaluate the right-hand side 

    subroutine sweet_f_eval(this, y, t, level_index, f, piece)
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

        if (piece == 1) then
            call ceval(y_sd_ptr%c_sweet_data_ptr, &
                    t,                         & 
                    this%ctx,                  &  
                    f_sd_ptr%c_sweet_data_ptr)
        else
            stop 'Bad value for piece in sweet_f_eval'
        end if 

    end subroutine sweet_f_eval

  
    ! solve for y in the implicit system + update f

    subroutine sweet_f_comp(this, y, t, dtq, rhs, level_index, f, piece)
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

        if (piece == 2) then
            call ccomp(y_sd_ptr%c_sweet_data_ptr,   & 
                    t,                           & 
                    dtq,                         & 
                    rhs_sd_ptr%c_sweet_data_ptr, &
                    this%ctx,                    & 
                    f_sd_ptr%c_sweet_data_ptr)
        else
            stop 'Bad value for piece in sweet_f_comp'
        end if
            
    end subroutine sweet_f_comp

    subroutine sweet_sweeper_initialize(this, pf, level_index)
        class(sweet_sweeper_t), intent(inout)        :: this
        type(pf_pfasst_t),      intent(inout),target :: pf
        integer,                intent(in)           :: level_index
        
        ! call superclass initialize
        call this%imex_initialize(pf, level_index)
        
        this%implicit=.TRUE.
        this%explicit=.TRUE.
    end subroutine sweet_sweeper_initialize
  
    ! destructor

    subroutine sweet_sweeper_destroy(this, pf, level_index)
        class(sweet_sweeper_t), intent(inout) :: this
        type(pf_pfasst_t),   intent(inout),target :: pf
        integer,             intent(in)    :: level_index

        ! this is copy-pasted from LibPFASST, unsure about this
        type(pf_level_t), pointer  :: lev       !  Current level
        lev => pf%levels(level_index)           !  Assign level pointer

        ! need the following line since the "final" keyword is not supported by some (older) compilers
        ! it forces Fortran to destroy the parent class data structures
        call this%imex_destroy(pf, level_index)
    end subroutine sweet_sweeper_destroy

    subroutine sweet_sweeper_compute_dt(this, pf, level_index, t0, dt, flags)
        class(sweet_sweeper_t),         intent(inout) :: this
        type(pf_pfasst_t), target,      intent(inout) :: pf
        integer,                        intent(in   ) :: level_index
        real(pfdp),                     intent(in   ) :: t0
        real(pfdp),                     intent(inout) :: dt
        integer, optional,              intent(in   ) :: flags

        type(pf_level_t),    pointer :: lev
        lev => pf%levels(level_index)   !!  Assign level pointer
        !  Do nothing now (copy-pasted from pf_imex_sweeper)
        return
    end subroutine sweet_sweeper_compute_dt

end module feval_module

