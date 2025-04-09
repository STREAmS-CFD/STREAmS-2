module catalyst2_conduit_wrappers

    use, intrinsic :: iso_fortran_env
    use iso_c_binding
    use catalyst_api
    use catalyst_conduit

#ifdef SINGLE_PRECISION
    integer, parameter, private :: rkind = REAL32
#else
    integer, parameter, private :: rkind = REAL64
#endif

contains

    subroutine catalyst_conduit_node_set_path_float_32_64(cnode, path, val) 
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        character(*), intent(IN) :: path
        real(rkind), value, intent(IN) :: val

#ifdef SINGLE_PRECISION
        call catalyst_conduit_node_set_path_float32(cnode, path, val)
#else
        call catalyst_conduit_node_set_path_float64(cnode, path, val)
#endif
    endsubroutine catalyst_conduit_node_set_path_float_32_64

    subroutine catalyst_conduit_node_set_external_float_32_64_ptr(cnode, data, num_elements) 
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        real(rkind), intent (IN), dimension (*) :: data
        integer(C_SIZE_T), value, intent(in) :: num_elements

#ifdef SINGLE_PRECISION
        call catalyst_conduit_node_set_external_float32_ptr(cnode, data, num_elements) 
#else
        call catalyst_conduit_node_set_external_float64_ptr(cnode, data, num_elements)
#endif
    endsubroutine catalyst_conduit_node_set_external_float_32_64_ptr

endmodule catalyst2_conduit_wrappers
