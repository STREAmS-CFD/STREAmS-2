program main

    use reader
    use postpro_cha
    use postpro_bl
    use postpro_chacurv
    use postpro_airfoil
    use postpro_ramp
    implicit none
    logical :: dir_exist

    call execute_command_line('mkdir -p POSTPRO/')

    call read_input
    call read_stat

    if (flow_init==0) call stats1d !channel flow statistics

    if (flow_init==1 .or. flow_init==2) then
        if(grid_dim == 1) then
            call stats2d !boundary layer flow statistics
        endif
        if(grid_dim == 2) then
            call stats2d_ramp !curved boundary layer/ramp flow statistics
        endif
    endif

    if (flow_init==4) call stats1d_curv ! curved channel

    if (flow_init==5) then  ! airfoil
        call forces_airfoil 
        call stats2d_airfoil
    endif

    if(save_plot3d > 0) call save_plot3d_fields()
end program main
