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
     
     subroutine c_sweet_data_restrict(io_y_coarse, i_y_fine, i_level_coarse, i_level_fine, i_t) & 
          bind(c, name="c_sweet_data_restrict")
       use iso_c_binding
       type(c_ptr),      value :: io_y_coarse, i_y_fine
       integer,          value :: i_level_coarse, i_level_fine
       double precision, value :: i_t
     end subroutine c_sweet_data_restrict

     subroutine c_sweet_data_interpolate(io_y_fine, i_y_coarse, i_level_fine, i_level_coarse, i_t) &
          bind(c, name="c_sweet_data_interpolate")
       use iso_c_binding
       type(c_ptr),      value :: io_y_fine, i_y_coarse
       integer,          value :: i_level_fine, i_level_coarse
       double precision, value :: i_t
     end subroutine c_sweet_data_interpolate
     
  end interface
  
contains

  ! interpolation function

  subroutine sweet_data_interpolate(this, levelF, levelG, qF, qG, t)
    class(sweet_level_t), intent(inout) :: this
    class(pf_level_t),       intent(inout) :: levelF, levelG
    class(pf_encap_t),       intent(inout) :: qF, qG
    real(pfdp),              intent(in)    :: t
    
    select type(qF)
    type is (sweet_data_encap_t)
       select type(qG)
       type is (sweet_data_encap_t)

          call c_sweet_data_interpolate(qF%c_sweet_data_ptr, & 
                                        qG%c_sweet_data_ptr, & 
                                        levelF%level,        &
                                        levelG%level,         &
                                        t)

       class default
          stop "TYPE ERROR"
       end select
    class default
       stop "TYPE ERROR"
    end select

  end subroutine sweet_data_interpolate

  
  ! restriction function

  subroutine sweet_data_restrict(this, levelF, levelG, qF, qG, t)
    class(sweet_level_t), intent(inout) :: this
    class(pf_level_t),       intent(inout) :: levelF, levelG
    class(pf_encap_t),       intent(inout) :: qF, qG
    real(pfdp),              intent(in)    :: t
    
    select type(qF)
    type is (sweet_data_encap_t)
       select type(qG)
       type is (sweet_data_encap_t)
          
          call c_sweet_data_restrict(qG%c_sweet_data_ptr, &
                                     qF%c_sweet_data_ptr, & 
                                     levelG%level,        & 
                                     levelF%level,        &
                                     t)

       class default
          stop "TYPE ERROR"
       end select
    class default
       stop "TYPE ERROR"
    end select

  end subroutine sweet_data_restrict

end module transfer_module

