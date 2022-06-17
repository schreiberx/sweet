module main_module
  use iso_c_binding 
  use encap_module
  use transfer_module
  use feval_module
  use hooks_module
  use pf_mod_rkstepper
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
    enddo
    ! LibPFASST supports more quadrature types but we do not support them here
    if (tmp(1:nl)      .eq. 'SDC_GAUSS_LOBATTO') then
       qtype = SDC_GAUSS_LOBATTO
    else if (tmp(1:nl) .eq. 'SDC_GAUSS_LEGENDRE') then
       qtype = SDC_GAUSS_LEGENDRE
    else
       qtype = -1
    endif
  end function translate_qtype

  ! main Fortran routine calling LibPFASST

  subroutine fmain(user_ctx_ptr,                                                           & ! user-defined context
                   nlevs, niters, nsweeps_coarse, nnodes, qtype_name, qnl, use_rk_stepper, & ! LibPFASST parameters
                   nfields, nvars_per_field,                                               & ! SWEET parameters
                   t_max, dt                                                               & ! timestepping parameters
                   ) bind (c, name='fmain')
    use mpi

    type(c_ptr),                 value       :: user_ctx_ptr
    integer                                  :: nlevs, niters, nsweeps_coarse, nnodes(nlevs), nvars(nlevs), shape(nlevs),   &
                                                nfields, nvars_per_field(nlevs), nsteps, level, qnl, qtype, use_rk_stepper, &
                                                ierror, num_procs, my_id, mpi_stat
    logical                                  :: use_no_left_q
    character(c_char)                        :: qtype_name
    real(c_double)                           :: t, t_max, dt
    class(pf_factory_t),         allocatable :: factory
    class(sweet_data_factory_t), pointer     :: sd_factory_ptr
    class(sweet_sweeper_t),      pointer     :: sweet_sweeper_ptr
    class(sweet_stepper_t),      pointer     :: sweet_stepper_ptr
    type(pf_comm_t)                          :: pf_comm
    type(pf_pfasst_t)                        :: pf

    real(c_double),             allocatable  :: z(:), y(:)
    real(c_double)                           :: val

    ! create the mpi and pfasst objects
     call pf_mpi_create(pf_comm, MPI_COMM_WORLD);
     call pf_pfasst_create(pf, pf_comm, nlevels=nlevs)

     call mpi_comm_rank(MPI_COMM_WORLD, my_id, ierror)
     call mpi_comm_size(MPI_COMM_WORLD, num_procs, ierror)

     if (my_id == 0) then
        print *, 'Number of Processors: ', num_procs
     end if
     
     ! timestepping parameters
     t      = 0
     nsteps = int(t_max/dt)

     ! LibPFASST parameters
     pf%nlevels           = nlevs                         ! number of SDC levels
     pf%niters            = niters                        ! number of SDC iterations
     if (use_rk_stepper == 1) then
        pf%use_rk_stepper = .true.                        ! replace the SDC sweeps with RK steps
        if (my_id == 0) then
           print *, "++++++++++++ use_rk_stepper == true  +++++++++++" 
        end if
     else
        pf%use_rk_stepper = .false.                       ! use the SDC sweeps
        if (my_id == 0) then
           print *, "++++++++++++ use_rk_stepper == false +++++++++++"
        end if
     end if
     pf%save_timings      = 1                             ! output the timings in fort.601 file
     qtype                = translate_qtype(qtype_name, & ! select the type of nodes
                                            qnl)

     !pf%abs_res_tol = 0.00000005 

     if (nlevs == 3) then
        nvars = [nfields*nvars_per_field(1), &
                 nfields*nvars_per_field(2), &
                 nfields*nvars_per_field(3)]    ! number of degrees of freedom for the levels
     else if (nlevs == 2) then
        nvars = [nfields*nvars_per_field(1), &
                 nfields*nvars_per_field(2)]    ! number of degrees of freedom for the levels
     else if (nlevs == 1) then
        nvars = [nfields*nvars_per_field(1)]    ! number of degrees of freedom for the levels
     else 
        stop 'This number of levels is not supported yet'
     end if

     ! loop over levels to initialize level-specific data structures
     do level = 1, pf%nlevels
       ! define the level id
       pf%levels(level)%index = level

       ! use the number of degrees of freedom as both the lev_shape and the buffer length of each level
       ! from what I understand, lev_shape refers to a collection of general variables relevant to the level
       ! buffer length: the amount of doubles to send in an MPI call
       call pf_level_set_size(pf, level, [nvars(level)], nvars(level))

       ! define the number of internal rk time steps
       pf%nsteps_rk = 1

       ! define number of sweeps for each level
       if (pf%nlevels == 1) then
          
          pf%levels(level)%nsweeps      = 1
          pf%levels(level)%nsweeps_pred = 1
             
       else

          if (level > 1) then
             
             pf%levels(level)%nsweeps      = 1             
             pf%levels(level)%nsweeps_pred = 1

          else
             pf%levels(level)%nsweeps      = nsweeps_coarse

             if (num_procs .eq. 1) then 
                pf%levels(level)%nsweeps_pred = 0 
             else
                pf%levels(level)%nsweeps_pred = 1
             end if
          end if
       end if

       ! allocate space for the levels

       !allocate(pf%levels%lev_shape(level))
       pf%levels(level)%lev_shape = nvars(level)
       
       ! define the properties (number of degrees of freedom and number of SDC nodes)
       ! pf%levels(level)%nvars   = nvars(level)
       pf%levels(level)%nnodes  = nnodes(level)
       pf%levels(level)%Finterp = .false.

       ! allocate space for the objects at this level
       allocate(sweet_level_t::pf%levels(level)%ulevel)
       allocate(sweet_data_factory_t::pf%levels(level)%ulevel%factory)
       allocate(sweet_sweeper_t::pf%levels(level)%ulevel%sweeper)
       allocate(sweet_stepper_t::pf%levels(level)%ulevel%stepper)
       
       ! select the order of the stepper
       if (level == pf%nlevels) then
          pf%levels(level)%ulevel%stepper%order = 4
       else
          pf%levels(level)%ulevel%stepper%order = 2
       end if
       
       ! check the number of nodes for RK stepper
       if (pf%use_rk_stepper .eqv. .true.) then
          if (pf%levels(level)%ulevel%stepper%order      == 5 .and. nnodes(level) /= 9) then
             stop "invalid number of nodes for the RK stepper"
          else if (pf%levels(level)%ulevel%stepper%order == 4 .and. nnodes(level) /= 7) then
             stop "invalid number of nodes for the RK stepper"
          else if (pf%levels(level)%ulevel%stepper%order == 3 .and. nnodes(level) /= 5) then
             stop "invalid number of nodes for the RK stepper"
          else if (pf%levels(level)%ulevel%stepper%order == 2 .and. nnodes(level) /= 4) then
             stop "invalid number of nodes for the RK stepper"
          end if
       end if

       ! cast the object into sweet data objects
       sd_factory_ptr    => as_sweet_data_factory(pf%levels(level)%ulevel%factory)
       sweet_sweeper_ptr => as_sweet_sweeper(pf%levels(level)%ulevel%sweeper)    
       sweet_stepper_ptr => as_sweet_stepper(pf%levels(level)%ulevel%stepper)

       ! pass the pointer to sweet data context to LibPFASST
       sd_factory_ptr%ctx    = user_ctx_ptr
       sweet_sweeper_ptr%ctx = user_ctx_ptr
       sweet_stepper_ptr%ctx = user_ctx_ptr

       ! initialize the sweeper data
       sweet_sweeper_ptr%level           = level
       sweet_sweeper_ptr%nnodes          = nnodes(level)
       sweet_sweeper_ptr%sweep_niter     = 0
       sweet_sweeper_ptr%sweep_niter_max = pf%niters
       sweet_sweeper_ptr%dt              = dt
       
       sweet_stepper_ptr%level           = level
       sweet_stepper_ptr%nnodes          = nnodes(level)
       sweet_stepper_ptr%sweep_niter     = 0
       sweet_stepper_ptr%sweep_niter_max = pf%niters
       sweet_stepper_ptr%dt              = dt

    end do

    ! initialize the pfasst objects
    call pf_pfasst_setup(pf)

    !! initialize the state vector at each level
    
    ! first fine level
    call finitial(pf%levels(pf%nlevels)%ulevel%sweeper, & 
                     pf%levels(pf%nlevels)%Q(1),        & 
                     pf%levels(pf%nlevels)%q0,          &
                     t_max,                             &
                     dt)
    ! then coarser levels
    if (pf%nlevels > 1) then
       do level = pf%nlevels-1, 1, -1
       
          call pf%levels(level)%ulevel%restrict(pf%levels(level+1),      &
                                                pf%levels(level),        &
                                                pf%levels(level+1)%Q(1), &
                                                pf%levels(level)%Q(1),   &
                                                t)

          call pf%levels(level)%q0%copy(pf%levels(level)%Q(1))
       end do
    end if
    

     call pf_add_hook(pf,                &
                      -1,        &
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

     call pf_add_hook(pf,                 &
                      pf%nlevels,         &
                      PF_POST_BLOCK,       &
                      fecho_output_invariants)
    
    ! advance in time with libpfasst
    level = nlevs

    call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)

    call pf_pfasst_run(pf,                    & 
                       pf%levels(level)%Q(1), &
                       dt,                    &
                       t,                     &
                       nsteps)
    t = t + dt*nsteps

    call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)
    
    ! finalize the simulation (does nothing right now)
    do level = 1, pf%nlevels

       call ffinal(pf%levels(level)%ulevel%sweeper,   & 
                   pf%levels(level)%Q(nnodes(level)), &
                   nnodes(level),                     &
                   pf%niters)

    end do

    call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)

    ! release memory
    call pf_pfasst_destroy(pf)

    !print *, 'all the memory was released'

  end subroutine fmain

end module main_module

