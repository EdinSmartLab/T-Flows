!==============================================================================!
  subroutine Control_Mod_Switch_To_Root()
!------------------------------------------------------------------------------!
!   Switch the control file to root.                                           !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  control_file_unit = root_control_file_unit

  end subroutine
