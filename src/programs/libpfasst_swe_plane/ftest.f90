module test_module
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

  public :: ftest

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

  subroutine ftest(user_ctx_ptr, num_levs, nvars) bind (c, name='ftest')
    use mpi

    type(c_ptr),                 value       :: user_ctx_ptr
    integer                                  :: num_levs, nnodes, level, nvars, shape(1)

    real(c_double)                           :: nnorm_s, nnorm_v, t, dt, t1, t2

    real(pfdp),                  allocatable :: z(:)

    class(pf_factory_t),         allocatable :: factory
    class(pf_encap_t),           allocatable :: encap_s_1
    class(pf_encap_t),           allocatable :: encap_s_2
    class(pf_encap_t),           allocatable :: encap_v_1(:)
    class(pf_encap_t),           allocatable :: encap_v_2(:)

    class(sweet_data_factory_t), pointer :: sd_factory_ptr
    class(sweet_data_encap_t),   pointer :: sd_encap_ptr   
    
    level  = 1
    nnodes = 2
    
    t1 = 1.d0

    ! allocate factory 

    allocate(sweet_data_factory_t::factory)
    sd_factory_ptr      => as_sweet_data_factory(factory)
    sd_factory_ptr%ctx  = user_ctx_ptr

    ! test create factory functions

    call sd_factory_ptr%create_single(encap_s_1, level, SDC_KIND_SOL_FEVAL, 3*nvars, shape)
    call sd_factory_ptr%create_single(encap_s_2, level, SDC_KIND_SOL_FEVAL, 3*nvars, shape)
    call sd_factory_ptr%create_array(encap_v_1, nnodes, level, SDC_KIND_SOL_FEVAL, 3*nvars, shape)
    call sd_factory_ptr%create_array(encap_v_2, nnodes, level, SDC_KIND_SOL_FEVAL, 3*nvars, shape)

    ! ! test setval, norm, and copy operators

    call encap_s_1%setval(t1)
    call encap_s_2%copy(encap_s_1)
    call encap_v_1(1)%setval(t1)
    call encap_v_2(1)%copy(encap_v_1(1))
    call encap_v_1(2)%setval(t1)
    call encap_v_2(2)%copy(encap_s_2)


    nnorm_s = encap_s_2%norm()
    nnorm_v = encap_v_2(1)%norm()
    
    print *, abs(nnorm_s - t1)
    print *, abs(nnorm_v - t1)
    
    if (abs(nnorm_s - t1) .gt. 1.d-20) then
       print *,'test1 norm is bad, ABS(norm - t1) =',ABS(nnorm_s-t1)
    endif
    
    if (abs(nnorm_v - t1) .gt. 1.d-20) then
       print *,'test1 norm is bad, ABS(norm - t1) =',ABS(nnorm_v-t1)
    endif

    ! test axpy operator
    
    t2 = 2.5d0
    call encap_s_2%axpy(t2,encap_s_1)
    call encap_v_2(1)%axpy(t2,encap_v_1(1))
    nnorm_s = encap_s_2%norm()
    nnorm_v = encap_v_2(1)%norm()

    if (abs(nnorm_s - t1*(1.d0+t2)) .gt. 1.d-20) then
       print *,'test2 norm is bad, ABS(norm - t1*(1+t2)) =',ABS(nnorm_s - t1*(1.d0+t2))
    endif

    if (abs(nnorm_v - t1*(1.d0+t2)) .gt. 1.d-20) then
       print *,'test2 norm is bad, ABS(norm - t1*(1+t2)) =',ABS(nnorm_v - t1*(1.d0+t2))
    endif
    
    print *, abs(nnorm_s - t1*(1.d0+t2))
    print *, abs(nnorm_v - t1*(1.d0+t2))

    ! test pack and unpack operators

    call encap_s_2%axpy(t2,encap_s_1)
    call encap_v_2(2)%axpy(t2,encap_v_1(1))

    allocate(z(3*nvars))
    call encap_s_2%pack(z)
    call encap_s_1%unpack(z)

    if (abs(encap_s_2%norm() - encap_s_1%norm()) .gt. 1.d-20) then
       print *,'test3 norm is bad, ABS(encap_s_2%norm() - encap_s_1%norm()) =',ABS(encap_s_2%norm() - encap_s_1%norm())
    endif

    print *, abs(encap_s_2%norm() - encap_s_1%norm())

    call encap_v_2(2)%pack(z)
    call encap_v_2(1)%unpack(z)

    if (abs(encap_v_2(2)%norm() - encap_v_2(1)%norm()) .gt. 1.d-20) then
       print *,'test4 norm is bad, ABS(encap_v_2(2)%norm() - encap_v_2(1)%norm()) =',ABS(encap_v_2(2)%norm() - encap_v_2(1)%norm())
    endif

    ! test destructors

    call factory%destroy_single(encap_s_1, level, SDC_KIND_SOL_FEVAL, 3*nvars, shape)
    call factory%destroy_single(encap_s_2, level, SDC_KIND_SOL_FEVAL, 3*nvars, shape)

    call factory%destroy_array(encap_v_1, nnodes, level, SDC_KIND_SOL_FEVAL, 3*nvars, shape)
    call factory%destroy_array(encap_v_2, nnodes, level, SDC_KIND_SOL_FEVAL, 3*nvars, shape)

    deallocate(factory)

  end subroutine ftest

end module test_module

