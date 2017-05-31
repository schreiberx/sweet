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
     
     subroutine c_sweet_data_restrict(Y_G, Y_F, t) bind(c, name="c_sweet_data_restrict")
       use iso_c_binding
       type(c_ptr),      value :: Y_G, Y_F
       double precision, value :: t
     end subroutine c_sweet_data_restrict


     subroutine c_sweet_data_interpolate(Y_F, Y_G, t) bind(c, name="c_sweet_data_interpolate")
       use iso_c_binding
       type(c_ptr),      value :: Y_F, Y_G
       double precision, value :: t
     end subroutine c_sweet_data_interpolate
     
  end interface
  
contains

  ! interpolation function

  subroutine sweet_data_interpolate(this, levelF, levelG, qF, qG, t)
    class(sweet_level_t), intent(inout) :: this
    class(pf_level_t),  intent(inout)   :: levelF, levelG
    class(pf_encap_t),  intent(inout)   :: qF, qG
    real(pfdp),         intent(in)      :: t

    select type(qF)
    type is (sweet_data_encap_t)
       select type(qG)
       type is (sweet_data_encap_t)

          call c_sweet_data_interpolate(qF%c_sweet_data_ptr, & 
                                        qG%c_sweet_data_ptr, & 
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
    class(pf_level_t),  intent(inout)   :: levelF, levelG
    class(pf_encap_t),  intent(inout)   :: qF, qG
    real(pfdp),         intent(in)      :: t
    
    select type(qF)
    type is (sweet_data_encap_t)
       select type(qG)
       type is (sweet_data_encap_t)
          
          call c_sweet_data_restrict(qG%c_sweet_data_ptr, &
                                     qF%c_sweet_data_ptr, & 
                                     t)

       class default
          stop "TYPE ERROR"
       end select
    class default
       stop "TYPE ERROR"
    end select

  end subroutine sweet_data_restrict

end module transfer_module

