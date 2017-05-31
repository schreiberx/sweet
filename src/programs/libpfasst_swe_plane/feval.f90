module feval_module
  use iso_c_binding 
  use pf_mod_dtype
  use encap_module
  use pf_mod_imexQ
  implicit none
  
  type, extends(pf_imexQ_t) :: sweet_sweeper_t
     type(c_ptr) :: ctx = c_null_ptr ! c pointer to PlaneDataCtx/SphereDataCtx
   contains 
     procedure :: f_eval  => sweet_f_eval 
     procedure :: f_comp  => sweet_f_comp
     procedure :: destroy => sweet_sweeper_destroy
  end type sweet_sweeper_t
  
  ! prototypes of the C functions

  interface 
     subroutine cinitial(i_ctx, i_t, i_dt, o_Y) bind(c, name="cinitial")
       use iso_c_binding
       type(c_ptr),    value :: i_ctx, o_Y
       real(c_double), value :: i_t, i_dt 
     end subroutine cinitial

     subroutine cfinal(i_ctx, i_Y) bind(c, name="cfinal")
       use iso_c_binding
       type(c_ptr), value :: i_ctx, i_Y
     end subroutine cfinal

     subroutine creference(i_ctx, t, o_Y) bind(c, name="creference")
       use iso_c_binding
       real(c_double), value :: t
       type(c_ptr),    value :: i_ctx, o_Y
     end subroutine creference

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

     subroutine ceval_f3(i_Y, i_t, i_ctx, o_F3) bind(c, name="ceval_f3")
       use iso_c_binding
       type(c_ptr),    value :: i_Y, i_ctx, o_F3 
       real(c_double), value :: i_t
     end subroutine ceval_f3

     subroutine ccomp_f3(i_Y, i_t, i_dt, i_Rhs, i_ctx, o_F3) bind(c, name="ccomp_f3")
       use iso_c_binding
       type(c_ptr),    value :: i_Y, i_Rhs, i_ctx, o_F3
       real(c_double), value :: i_t, i_dt
     end subroutine ccomp_f3
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
                  t, &
                  dt, &
                  sd_ptr%c_sweet_data_ptr)

  end subroutine finitial


  subroutine ffinal(sweeper, sd)
    class(pf_sweeper_t),       intent(in)    :: sweeper
    class(pf_encap_t),         intent(inout) :: sd

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: sd_ptr

    sweet_sweeper_ptr => as_sweet_sweeper(sweeper)
    sd_ptr            => as_sweet_data_encap(sd)

    call cfinal(sweet_sweeper_ptr%ctx, & 
                sd_ptr%c_sweet_data_ptr)

  end subroutine ffinal

  
  subroutine freference(sweeper, t, sd)
    class(pf_sweeper_t),       intent(in)    :: sweeper
    real(pfdp),                intent(in)    :: t
    class(pf_encap_t),         intent(inout) :: sd

    class(sweet_sweeper_t),    pointer       :: sweet_sweeper_ptr
    class(sweet_data_encap_t), pointer       :: sd_ptr

    sweet_sweeper_ptr => as_sweet_sweeper(sweeper)
    sd_ptr            => as_sweet_data_encap(sd)

    call creference(sweet_sweeper_ptr%ctx, &
                    t, &
                    sd_ptr%c_sweet_data_ptr)

  end subroutine freference
  
  ! evaluate the right-hand side 

  subroutine sweet_f_eval(this, y, t, level, f, piece)
    class(sweet_sweeper_t),    intent(inout) :: this
    class(pf_encap_t),         intent(in)    :: y
    real(pfdp),                intent(in)    :: t
    integer,                   intent(in)    :: level
    class(pf_encap_t),         intent(inout) :: f
    integer,                   intent(in)    :: piece
    class(sweet_data_encap_t), pointer       :: y_sd_ptr
    class(sweet_data_encap_t), pointer       :: f_sd_ptr
    
    y_sd_ptr  => as_sweet_data_encap(y)
    f_sd_ptr  => as_sweet_data_encap(f)

    select case (piece)

       case (1) ! explicit rhs
          call ceval_f1(y_sd_ptr%c_sweet_data_ptr, &
                        t, & 
                        this%ctx, & 
                        f_sd_ptr%c_sweet_data_ptr)

       case (2) ! first implicit rhs
          call ceval_f2(y_sd_ptr%c_sweet_data_ptr, & 
                        t, & 
                        this%ctx, & 
                        f_sd_ptr%c_sweet_data_ptr)

       ! case (3) ! second implicit rhs
       !    call ceval_f3(y_sd_ptr%c_sweet_data_ptr, & 
       !                  t, & 
       !                  this%ctx, &
       !                  f_sd_ptr%c_sweet_data_ptr)

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
          call ccomp_f2(y_sd_ptr%c_sweet_data_ptr, & 
                        t, & 
                        dt, & 
                        rhs_sd_ptr%c_sweet_data_ptr, &
                        this%ctx, & 
                        f_sd_ptr%c_sweet_data_ptr)

    !    case (3) ! second implicit solve
    !       call ccomp_f3(y_sd_ptr%c_sweet_data_ptr, &
    !                     t, & 
    !                     dt, & 
    !                     rhs_sd_ptr%c_sweet_data_ptr, & 
    !                     this%ctx, & 
    !                     f_sd_ptr%c_sweet_data_ptr)

       case DEFAULT
          print *, 'Piece argument in f_comp can only be 2'! or 3'
          call exit(0)
    end select
         
  end subroutine sweet_f_comp

  
  ! destructor

  subroutine sweet_sweeper_destroy(this, lev)
    class(sweet_sweeper_t), intent(inout) :: this
    class(pf_level_t),      intent(inout) :: lev

    ! need the follwing line since the "final" keyword is not supported by some (older) compilers
    ! it forces Fortran to destroy the parent class data structures
    call this%imexQ_destroy(lev) 

  end subroutine sweet_sweeper_destroy

end module feval_module

