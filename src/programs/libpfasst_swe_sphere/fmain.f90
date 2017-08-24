module main_module
  use iso_c_binding 
  use encap_module
!  use transfer_module
!  use feval_module
!  use hooks_module
  use pf_mod_misdc
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
    do i=1,nl
       tmp(i:i) = name(i)
    enddo
    if (tmp(1:nl)      .eq. 'SDC_GAUSS_LOBATTO') then
       qtype = SDC_GAUSS_LOBATTO
    else if (tmp(1:nl) .eq. 'SDC_GAUSS_RADAU') then
       qtype = SDC_GAUSS_RADAU
    else if (tmp(1:nl) .eq. 'SDC_CLENSHAW_CURTIS') then
       qtype = SDC_CLENSHAW_CURTIS
    else if (tmp(1:nl) .eq. 'SDC_UNIFORM') then
       qtype = SDC_UNIFORM
    else if (tmp(1:nl) .eq. 'SDC_GAUSS_LEGENDRE') then
       qtype = SDC_GAUSS_LEGENDRE
    else if (tmp(1:nl) .eq. 'SDC_PROPER_NODES') then
       qtype = SDC_PROPER_NODES
    else if (tmp(1:nl) .eq. 'SDC_COMPOSITE_NODES') then
       qtype = SDC_COMPOSITE_NODES
    else if (tmp(1:nl)      .eq. 'SDC_GAUSS_LOBATTO_NL') then
       qtype = SDC_GAUSS_LOBATTO + SDC_NO_LEFT
    else if (tmp(1:nl) .eq. 'SDC_GAUSS_RADAU_NL') then
       qtype = SDC_GAUSS_RADAU + SDC_NO_LEFT
    else if (tmp(1:nl) .eq. 'SDC_CLENSHAW_CURTIS_NL') then
       qtype = SDC_CLENSHAW_CURTIS + SDC_NO_LEFT
    else if (tmp(1:nl) .eq. 'SDC_UNIFORM_NL') then
       qtype = SDC_UNIFORM + SDC_NO_LEFT
    else if (tmp(1:nl) .eq. 'SDC_GAUSS_LEGENDRE_NL') then
       qtype = SDC_GAUSS_LEGENDRE + SDC_NO_LEFT
    else if (tmp(1:nl) .eq. 'SDC_PROPER_NODES_NL') then
       qtype = SDC_PROPER_NODES + SDC_NO_LEFT
    else if (tmp(1:nl) .eq. 'SDC_COMPOSITE_NODES_NL') then
       qtype = SDC_COMPOSITE_NODES + SDC_NO_LEFT
    else
       qtype = -1
    endif
  end function translate_qtype

  ! main Fortran routine calling LibPFASST

  subroutine fmain(                                        &
                   user_ctx_ptr,                           & ! user-defined context
                   nlevs, niters, nnodes, qtype_name, qnl, & ! LibPFASST parameters
                   nfields, nvars_per_field,               & ! SWEET parameters
                   t_max, dt                               & ! timestepping parameters
                   ) bind (c, name='fmain')
    use mpi

    type(c_ptr),                 value       :: user_ctx_ptr
    integer                                  :: nlevs, niters, nnodes(nlevs), nvars(nlevs), shape(nlevs), &
                                                nfields, nvars_per_field(nlevs), nsteps, level, qnl, qtype, &
                                                ierr, num_procs, my_id, mpi_stat
    character(c_char)                        :: qtype_name
    real(c_double)                           :: t, t_max, dt
    class(pf_factory_t),         allocatable :: factory
    class(sweet_data_factory_t), pointer     :: sd_factory_ptr
!    class(sweet_sweeper_t),      pointer     :: sweet_sweeper_ptr
    type(pf_comm_t)                          :: pf_comm
    type(pf_pfasst_t)                        :: pf

    real(c_double),             allocatable  :: z(:)

    ! create the mpi and pfasst objects
     call pf_mpi_create(pf_comm,    & 
                        MPI_COMM_WORLD);
     call pf_pfasst_create(pf,      & 
                           pf_comm, & 
                           nlevs)

     call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
     call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

 !     if (my_id == 0) then
!         print *, 'Number of Processors: ', num_procs
!      end if
!      print *, 'My Id = ', my_id

!      ! timestepping parameters
!      t      = 0
!      nsteps = int(t_max/dt)

!      ! LibPFASST parameters
!      pf%nlevels    = nlevs   ! number of SDC levels
!      pf%niters     = niters  ! number of SDC iterations
!      pf%pipeline_G = .true.  ! pipeline the coarse prediction sweeps
!      qtype         = translate_qtype(qtype_name, & ! select the type of nodes
!                                      qnl)

!      if (nlevs == 3) then
!         nvars = [nfields*nvars_per_field(1), &
!                  nfields*nvars_per_field(2), &
!                  nfields*nvars_per_field(3)]    ! number of degrees of freedom for the levels
!      else if (nlevs == 2) then
!         nvars = [nfields*nvars_per_field(1), &
!                  nfields*nvars_per_field(2)]    ! number of degrees of freedom for the levels
!      else if (nlevs == 1) then
!         nvars = [nfields*nvars_per_field(1)]    ! number of degrees of freedom for the levels
!      else 
!         stop 'This number of levels is not supported yet'
!      end if

!      ! loop over levels to initialize level-specific data structures
!      do level = 1, pf%nlevels

!        ! define the level id
!        pf%levels(level)%level = level ! confusing!

!        ! define number of sweeps for each level
!        pf%levels(level)%nsweeps_pred = 1
!        if (level > 1) then
!           pf%levels(level)%nsweeps   = 1
!        else
!           pf%levels(level)%nsweeps   = 2
!        end if

!        ! allocate space for the levels
!        allocate(pf%levels(level)%shape(level))
!        pf%levels(level)%shape(level) = nvars(level)
       
!        ! define the properties (number of degrees of freedom and number of SDC nodes)
!        pf%levels(level)%nvars  = nvars(level)
!        pf%levels(level)%nnodes = nnodes(level)
       
!        ! allocate space for the objects at this level
!        allocate(sweet_level_t::pf%levels(level)%ulevel)
!        allocate(sweet_data_factory_t::pf%levels(level)%ulevel%factory)
!        allocate(sweet_sweeper_t::pf%levels(level)%ulevel%sweeper)

!        ! cast the object into sweet data objects
!        sd_factory_ptr    => as_sweet_data_factory(pf%levels(level)%ulevel%factory)
!        sweet_sweeper_ptr => as_sweet_sweeper(pf%levels(level)%ulevel%sweeper)    

!        ! pass the pointer to sweet data context to LibPFASST
!        sd_factory_ptr%ctx    = user_ctx_ptr
!        sweet_sweeper_ptr%ctx = user_ctx_ptr
       
!        ! initialize the sweeper data
!        sweet_sweeper_ptr%level           = level
!        sweet_sweeper_ptr%nnodes          = nnodes(level)
!        sweet_sweeper_ptr%sweep_niter     = 0
!        sweet_sweeper_ptr%sweep_niter_max = pf%niters
!        sweet_sweeper_ptr%dt              = dt
       
!     end do

!     ! initialize the mpi and pfasst objects
!     call pf_mpi_setup(pf_comm, & 
!                       pf)
!     call pf_pfasst_setup(pf)
    
!     ! initialize the state vector at each level
!     do level = 1, pf%nlevels

!        call finitial(pf%levels(level)%ulevel%sweeper, & 
!                      pf%levels(level)%Q(1),           & 
!                      pf%levels(level)%q0,             &
!                      t,                               &
!                      dt)

!     end do

! !    allocate(z(nvars(pf%nlevels)))
! !    call pf%levels(pf%nlevels)%q0%unpack(z)
! !    call pf%levels(pf%nlevels)%q0%pack(z)

  
!    ! define the hooks to output data to the terminal (residual and error)
!     ! call pf_add_hook(pf,            &
!     !                  pf%nlevels,    &
!     !                  PF_POST_SWEEP, &
!     !                  fecho_error)
!     call pf_add_hook(pf,             &
!                      -1,             &
!                      PF_POST_SWEEP,  &
!                      fecho_residual)   
!     ! call pf_add_hook(pf,           &
!     !                  pf%nlevels,   &
!     !                  PF_POST_STEP, &
!     !                  fecho_interpolation_errors)
!     call pf_add_hook(pf,                 &
!                      pf%nlevels,         &
!                      PF_POST_ITERATION,  &
!                      fecho_output_solution)

!     ! advance in time with libpfasst
!     level = nlevs

!     call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)

!     call pf_pfasst_run(pf,                    & 
!                        pf%levels(level)%Q(1), &
!                        dt,                    &
!                        t,                     &
!                        nsteps)
!     t = t + dt*nsteps

!     ! finalize the simulation (does nothing right now)
!     do level = 1, pf%nlevels

!        call ffinal(pf%levels(level)%ulevel%sweeper,   & 
!                    pf%levels(level)%Q(nnodes(level)), &
!                    nnodes(level),                     &
!                    pf%niters)

!     end do

!     call MPI_BARRIER(MPI_COMM_WORLD,mpi_stat)

    ! release memory
    call pf_pfasst_destroy(pf)
    call pf_mpi_destroy(pf_comm)

    print *, 'all the memory was released'

  end subroutine fmain

end module main_module

