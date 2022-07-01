module main_module
    use iso_c_binding 
    use encap_module
    use feval_module
    use hooks_module
    use transfer_module
    use pf_mod_parallel
    use pf_mod_mpi
    use pfasst
    implicit none

    public :: fmain

contains
    function translate_qtype(name, nl) result (qtype)
        integer :: nl
        character(c_char), intent(in) :: name(nl)
        character*72 :: tmp
        integer :: qtype, i
        logical :: use_no_left_q
        do i=1,nl
            tmp(i:i) = name(i)
        end do
        ! LibPFASST supports more quadrature types but we do not support them here
        if (tmp(1:nl)      .eq. 'SDC_GAUSS_LOBATTO') then
            qtype = SDC_GAUSS_LOBATTO
        else if (tmp(1:nl) .eq. 'SDC_GAUSS_LEGENDRE') then
            qtype = SDC_GAUSS_LEGENDRE
        else
            qtype = -1
        end if
    end function translate_qtype

    ! main Fortran routine calling LibPFASST

    subroutine fmain(user_ctx_ptr,                                                           & ! user-defined context
                    niters, nsweeps, nnodes, qtype_name, qnl, use_rk_stepper, & ! LibPFASST parameters
                    nfields, nvars_per_field,                                               & ! SWEET parameters
                    t_max, dt                                                               & ! timestepping parameters
                    ) bind (c, name='fmain')
        use mpi

        type(c_ptr),                 value       :: user_ctx_ptr
        integer                                  :: niters, nsweeps, nnodes, nvars,   &
                                                    nfields, nvars_per_field, nsteps, qnl, qtype, use_rk_stepper, &
                                                    ierror, num_procs, my_id, mpi_stat
        logical                                  :: use_no_left_q
        character(c_char)                        :: qtype_name
        real(c_double)                           :: t, t_max, dt
        class(pf_factory_t),         allocatable :: factory
        class(sweet_data_factory_t), pointer     :: sd_factory_ptr
        class(sweet_sweeper_t),      pointer     :: sweet_sweeper_ptr
        type(pf_comm_t)                          :: pf_comm
        type(pf_pfasst_t)                        :: pf

        real(c_double),             allocatable  :: z(:), y(:)
        real(c_double)                           :: val

        print *, 'in fmain'

        ! create the mpi and pfasst objects
        call pf_mpi_create(pf_comm, MPI_COMM_WORLD);
        print *, 'created mpi object'
        call pf_pfasst_create(pf, pf_comm, nlevels=1)
        print *, 'created pfasst object'

        call mpi_comm_rank(MPI_COMM_WORLD, my_id, ierror)
        call mpi_comm_size(MPI_COMM_WORLD, num_procs, ierror)

        if (my_id == 0) then
            print *, 'Number of Processors: ', num_procs
        end if
        
        ! timestepping parameters
        t      = 0
        nsteps = int(t_max/dt)

        ! LibPFASST parameters
        pf%niters            = niters                        ! number of SDC iterations
        pf%save_timings      = 1                             ! output the timings in fort.601 file
        pf%qtype             = translate_qtype(qtype_name, & ! select the type of nodes
                                            qnl)

        nvars = nfields*nvars_per_field    ! number of degrees of freedom for the levels

        ! initialize level-specific data structures
        pf%levels(1)%index = 1
        
        call pf_level_set_size(pf, 1, [nvars], nvars)

        ! define the number of internal rk time steps
        pf%nsteps_rk = 1

        pf%levels(1)%nsweeps      = nsweeps
        pf%levels(1)%nsweeps_pred = 1
                

        ! allocate space for the levels
        pf%levels(1)%lev_shape = nvars
        
        ! define the properties (number of degrees of freedom and number of SDC nodes)
        pf%levels(1)%nnodes  = nnodes
        pf%levels(1)%Finterp = .false.

        ! allocate space for the objects at this level
        allocate(sweet_level_t::pf%levels(1)%ulevel)
        allocate(sweet_data_factory_t::pf%levels(1)%ulevel%factory)
        allocate(sweet_sweeper_t::pf%levels(1)%ulevel%sweeper)

        ! cast the object into sweet data objects
        sd_factory_ptr    => as_sweet_data_factory(pf%levels(1)%ulevel%factory)
        sweet_sweeper_ptr => as_sweet_sweeper(pf%levels(1)%ulevel%sweeper)    

        ! pass the pointer to sweet data context to LibPFASST
        sd_factory_ptr%ctx    = user_ctx_ptr
        sweet_sweeper_ptr%ctx = user_ctx_ptr

        ! initialize the sweeper data
        sweet_sweeper_ptr%level           = 1
        sweet_sweeper_ptr%nnodes          = nnodes
        sweet_sweeper_ptr%sweep_niter     = 0
        sweet_sweeper_ptr%sweep_niter_max = pf%niters
        sweet_sweeper_ptr%dt              = dt

        ! initialize the pfasst objects
        call pf_pfasst_setup(pf)

        !! initialize the state vector
        call finitial(pf%levels(pf%nlevels)%ulevel%sweeper, & 
                      pf%levels(pf%nlevels)%Q(1),           & 
                      pf%levels(pf%nlevels)%q0,             &
                      t_max,                                &
                      dt)
        

        call pf_add_hook(pf,                &
                         -1,                &
                         PF_POST_ITERATION, &
                         fecho_residual)
        if (num_procs .eq. 1) then
            call pf_add_hook(pf,                 &
                             -1,                 &
                             PF_POST_ITERATION,  &
                             fecho_output_solution)
        else
            call pf_add_hook(pf,                 &
                             -1,                 &
                             PF_POST_SWEEP,      &
                             fecho_output_solution)
        end if

        call pf_add_hook(pf,                  &
                         pf%nlevels,          &
                         PF_POST_BLOCK,       &
                         fecho_output_invariants)
        
        call pf_print_options(pf,un_opt=6)
        
        ! advance in time with libpfasst

        call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)

        call pf_pfasst_run(pf,                    & 
                           pf%levels(1)%Q(1), &
                           dt,                    &
                           t,                     &
                           nsteps)
        t = t + dt*nsteps

        call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)
        
        ! finalize the simulation (does nothing right now)
        call ffinal(pf%levels(1)%ulevel%sweeper,   & 
                    pf%levels(1)%Q(nnodes), &
                    nnodes,                     &
                    pf%niters)


        call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)

        ! release memory
        call pf_pfasst_destroy(pf)

    end subroutine fmain

end module main_module

