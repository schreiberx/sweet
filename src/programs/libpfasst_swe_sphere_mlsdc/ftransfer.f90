module transfer_module
  use pf_mod_dtype
  use encap_module
  use feval_module
  implicit none

  type, extends(pf_user_level_t) :: sweet_level_t
   contains 
     procedure :: restrict    => sweet_data_restrict 
     procedure :: interpolate => sweet_data_interpolate
  end type sweet_level_t

  interface 

     ! prototypes of the C functions
     
     subroutine c_sweet_data_restrict(io_y_coarse, i_y_fine, i_level_coarse, i_level_fine, i_ctx, i_t) & 
          bind(c, name="c_sweet_data_restrict")
       use iso_c_binding
       type(c_ptr),      value :: io_y_coarse, i_y_fine, i_ctx
       integer,          value :: i_level_coarse, i_level_fine
       double precision, value :: i_t
     end subroutine c_sweet_data_restrict

     subroutine c_sweet_data_interpolate(io_y_fine, i_y_coarse, i_level_fine, i_level_coarse, i_ctx, i_t) &
          bind(c, name="c_sweet_data_interpolate")
       use iso_c_binding
       type(c_ptr),      value :: io_y_fine, i_y_coarse, i_ctx
       integer,          value :: i_level_fine, i_level_coarse
       double precision, value :: i_t
     end subroutine c_sweet_data_interpolate
     
  end interface
  
contains

  ! interpolation function

  subroutine sweet_data_interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(sweet_level_t),    intent(inout) :: this
    class(pf_level_t),       intent(inout) :: f_lev, c_lev ! fine and coarse levels
    class(pf_encap_t),       intent(inout) :: f_vec, c_vec ! fine and coarse vectors
    real(pfdp),              intent(in)    :: t
    integer, optional,       intent(in)    :: flags

    class(sweet_sweeper_t),  pointer       :: sweet_sweeper_ptr

    sweet_sweeper_ptr => as_sweet_sweeper(this%sweeper)

    select type(f_vec)
    type is (sweet_data_encap_t)
       select type(c_vec)
       type is (sweet_data_encap_t)

          call c_sweet_data_interpolate(f_vec%c_sweet_data_ptr,   &
                                        c_vec%c_sweet_data_ptr,   &
                                        f_lev%index-1,        & ! conversion to c++ indexing
                                        c_lev%index-1,        & ! conversion to c++ indexing
                                        sweet_sweeper_ptr%ctx, &
                                        t)

       class default
          stop "TYPE ERROR"
       end select
    class default
       stop "TYPE ERROR"
    end select

  end subroutine sweet_data_interpolate

  
  ! restriction function

  subroutine sweet_data_restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(sweet_level_t),    intent(inout) :: this
    class(pf_level_t),       intent(inout) :: f_lev, c_lev ! fine and coarse levels
    class(pf_encap_t),       intent(inout) :: f_vec, c_vec ! fine and coarse vectors
    real(pfdp),              intent(in)    :: t
    integer, optional,       intent(in)    :: flags
    
    class(sweet_sweeper_t),  pointer       :: sweet_sweeper_ptr

    sweet_sweeper_ptr => as_sweet_sweeper(this%sweeper)

    select type(f_vec)
    type is (sweet_data_encap_t)
       select type(c_vec)
       type is (sweet_data_encap_t)
          
          call c_sweet_data_restrict(c_vec%c_sweet_data_ptr,   &
                                     f_vec%c_sweet_data_ptr,   &
                                     c_lev%index-1,        & ! conversion to c++ indexing
                                     f_lev%index-1,        & ! conversion to c++ indexing
                                     sweet_sweeper_ptr%ctx, &
                                     t)

       class default
          stop "TYPE ERROR"
       end select
    class default
       stop "TYPE ERROR"
    end select

  end subroutine sweet_data_restrict

end module transfer_module

