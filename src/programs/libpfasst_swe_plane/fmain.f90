module main_module
  use iso_c_binding 
  use encap_module
  use transfer_module
  use feval_module
  use hooks_module
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

  subroutine fmain(user_ctx_ptr, num_levs, nvars) bind (c, name='fmain')
    use mpi

    type(c_ptr),                  value      :: user_ctx_ptr
    integer                                  :: num_levs, nnodes(1), qnl, iters, steps, nsteps, level, nvars, shape(1), kind, qtype
    character(c_char)                        :: qtype_name

    real(c_double)                           :: t, dt

    class(pf_factory_t),         allocatable :: factory
    class(pf_encap_t),           allocatable :: q0
    
    class(sweet_data_factory_t), pointer     :: sd_factory_ptr
    class(sweet_sweeper_t),      pointer     :: sweet_sweeper_ptr
    
    type(pf_comm_t)   :: pf_comm
    type(pf_pfasst_t) :: pf


    ! create the mpi and pfasst objects
     call pf_mpi_create(pf_comm, & 
                        MPI_COMM_SELF);
     call pf_pfasst_create(pf, & 
                           pf_comm, & 
                           num_levs)
     pf_comm%nproc = 1

     ! parameters definition
     steps      = 1 
     nsteps     = steps*pf_comm%nproc
     
     t          = 0.0_pfdp
     dt         = 0.0001_pfdp ! small time step for now since the integrator is explicit
     
     pf%nlevels = num_levs
     pf%niters  = 20 ! number of iterations hard coded for now
     
     qtype_name = 'SDC_GAUSS_LOBATTO' ! type of nodes hard coded for now
     qtype      = translate_qtype(qtype_name, & 
                                  qnl)
     nnodes(1)  = 6 ! number of nodes hard coded for now
 
     ! loop over levels (currently only one) to initialize level-specific data structures
     do level = 1, pf%nlevels

       ! trivial zero-order predictor
       pf%levels(level)%nsweeps_pred = 1
       
       ! allocate space for the levels
       allocate(pf%levels(level)%shape(1))
       pf%levels(level)%shape(1) = pf%levels(level)%nvars

       ! define the properties
       pf%levels(level)%nvars  = nvars 
       pf%levels(level)%nnodes = nnodes(level)
       
       ! allocate space for the objects at this level
       allocate(sweet_level_t::pf%levels(level)%ulevel)
       allocate(sweet_data_factory_t::pf%levels(level)%ulevel%factory)
       allocate(sweet_sweeper_t::pf%levels(level)%ulevel%sweeper)

       ! cast the object into sweet data objects
       sd_factory_ptr    => as_sweet_data_factory(pf%levels(level)%ulevel%factory)
       sweet_sweeper_ptr => as_sweet_sweeper(pf%levels(level)%ulevel%sweeper)    
       
       ! pass the pointer to sweet data context to libpfasst
       sd_factory_ptr%ctx    = user_ctx_ptr
       sweet_sweeper_ptr%ctx = user_ctx_ptr
       
    end do

    ! initialize the mpi and pfasst objects
    call pf_mpi_setup(pf_comm, & 
                      pf)
    call pf_pfasst_setup(pf)
    
    ! initialize the state vector at each level
    level = 1
    call sd_factory_ptr%create_single(q0, & 
                                      level, & 
                                      SDC_KIND_SOL_FEVAL, & 
                                      nvars, &
                                      shape)

    do level = 1, pf%nlevels

       call finitial(pf%levels(level)%ulevel%sweeper, & 
                     pf%levels(level)%Q(1), & 
                     pf%levels(level)%q0, &
                     t, &
                     dt)

    end do
    
    ! define the hooks to output data to the terminal (residual and error)
    call pf_add_hook(pf, & 
                     pf%nlevels, & 
                     PF_POST_STEP, &
                     fecho_error)
    call pf_add_hook(pf, &
                     -1, &
                     PF_POST_SWEEP, &
                     fecho_residual)
    

    ! advance in time with libpfasst
    level = 1
    call pf_pfasst_run(pf, & 
                       pf%levels(level)%Q(1), &
                       dt, &
                       t, &
                       nsteps)
    t = t + dt*nsteps

    ! finalize the simulation (does nothing right now)
    do level = 1, pf%nlevels
       call ffinal(pf%levels(level)%ulevel%sweeper, & 
                   pf%levels(level)%Q(nnodes(level)))
    end do
    
    ! release memory
    level = 1
    call sd_factory_ptr%destroy_single(q0, &
                                       level, & 
                                       SDC_KIND_SOL_FEVAL, &
                                       nvars, &
                                       shape)

    call pf_pfasst_destroy(pf)
    call pf_mpi_destroy(pf_comm)

  end subroutine fmain

end module main_module

