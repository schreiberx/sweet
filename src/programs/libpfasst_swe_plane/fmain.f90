module main_module
  use iso_c_binding 
  use encap_module
  use transfer_module
  use feval_module
  use hooks_module
  use pf_mod_misdc
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

  subroutine fmain(userCtx_ptr,num_levs,nodes,qtype_name,qnl,iters,steps,t,dt)  bind(c, name='fmain')
    use parallel

    type(c_ptr),value :: userCtx_ptr
    integer           :: num_levs, nodes(*), qnl, iters, steps
    character(c_char) :: qtype_name(qnl)
    real(c_double)    :: t, dt
    type(pf_comm_t)   :: pf_comm
    type(pf_pfasst_t) :: pf
    type(pf_encap_t), target :: encap

    integer     :: level, shape(1), kind, nvars, myproc, nsteps
    real(c_double) :: norm

    type(c_ptr) ::  q0 = c_null_ptr
    type(c_ptr) ::  qN = c_null_ptr
    type(c_ptr) :: rhs = c_null_ptr
    character*72 :: fname

    nvars = 1 ! Unused ...
    call c_encap_encap_create(encap)
    call pf_mpi_create(pf_comm, MPI_COMM_SELF);
    call pf_pfasst_create(pf, pf_comm, num_levs, fname='pfasst.nml')

    pf%qtype  = translate_qtype(qtype_name,qnl)
    pf%niters = iters

    do level = 1, pf%nlevels

       pf%levels(level)%nsweeps_pred = 1

       pf%levels(level)%nvars  = nvars ! Not actually used, but must be set
       pf%levels(level)%nnodes = nodes(level)
       pf%levels(level)%ctx    = userCtx_ptr

       pf%levels(level)%interpolate => c_encap_interpolate
       pf%levels(level)%restrict    => c_encap_restrict
       pf%levels(level)%encap       => encap

       !call pf_explicit_create(pf%levels(level)%sweeper, feval_f1)
       !call pf_imex_create(pf%levels(level)%sweeper, feval_f1, feval_f2, fcomp_f2)
       call pf_misdc_create(pf%levels(level)%sweeper, feval_f1, feval_f2, fcomp_f2, feval_f3, fcomp_f3)
    end do

    call pf_mpi_setup(pf_comm, pf)
    call pf_pfasst_setup(pf)

    do level = 1, pf%nlevels
       call finitial(pf%levels(level)%q0)
    end do

    call pf_add_hook(pf, pf%nlevels, PF_POST_STEP, fecho_error)
    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)

    nsteps = steps*pf_comm%nproc
    call pf_pfasst_run(pf, q0, dt, tend=t, nsteps=nsteps)
    t = t + dt*nsteps

    do level = 1, pf%nlevels
       call ffinal(pf%levels(level)%Q(nodes(level)))
    end do

    call pf_pfasst_destroy(pf)

  end subroutine fmain

end module main_module
