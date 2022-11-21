program main

 use reader
 use postpro_cha
 use postpro_bl
 implicit none
 logical :: dir_exist
!
 call execute_command_line('mkdir -p POSTPRO/')
!
 call read_input
 call read_stat
 if (flow_init==0) call stats1d !channel flow statistics
 if (flow_init>=1) call stats2d !boundary layer flow statistics
!

end program main
