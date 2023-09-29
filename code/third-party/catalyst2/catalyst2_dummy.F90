module catalyst_conduit

contains

    subroutine catalyst_conduit_node_info(cnode,cdest) 
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        type(C_PTR), value, intent(IN) :: cdest
    end subroutine catalyst_conduit_node_info

    subroutine catalyst_conduit_node_print(cnode) 
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
    end subroutine catalyst_conduit_node_print

    function catalyst_conduit_node_create() result(cnode)
        use iso_c_binding
        implicit none
        type(C_PTR) :: cnode
        cnode = c_null_ptr
    end function catalyst_conduit_node_create

    subroutine catalyst_conduit_node_set_path_int64(cnode, path, val)
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        character(*), intent(IN) :: path
        integer(8), value, intent(IN) :: val
    end subroutine catalyst_conduit_node_set_path_int64

    subroutine catalyst_conduit_node_set_path_float64(cnode, path, val)
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        character(*), intent(IN) :: path
        real(8), value, intent(IN) :: val
    end subroutine catalyst_conduit_node_set_path_float64

    subroutine catalyst_conduit_node_set_path_char8_str(cnode, path, val)
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        character(*), intent(IN) :: path
        character(*), intent(IN) :: val
    end subroutine catalyst_conduit_node_set_path_char8_str

    subroutine c_catalyst_conduit_node_set_path_char8_str(cnode, path, val) 
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        character(kind=C_CHAR), intent(IN) :: path(*)
        character(kind=C_CHAR), intent(IN) :: val(*)
    end subroutine c_catalyst_conduit_node_set_path_char8_str

    subroutine catalyst_conduit_node_set_external_float64_ptr(cnode, data, num_elements) 
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        real(8), intent (IN), dimension (*) :: data
        integer(C_SIZE_T), value, intent(in) :: num_elements
    end subroutine catalyst_conduit_node_set_external_float64_ptr

    subroutine catalyst_conduit_node_set_path_external_node(cnode, path, cother)
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        character(*), intent(IN) :: path
        type(C_PTR), value, intent(IN) :: cother
    end subroutine catalyst_conduit_node_set_path_external_node

    subroutine catalyst_conduit_node_set_path_int32(cnode, path, val)
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        character(*), intent(IN) :: path
        integer(4), value, intent(IN) :: val
    end subroutine catalyst_conduit_node_set_path_int32

    subroutine catalyst_conduit_node_destroy(cnode) 
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
    end subroutine catalyst_conduit_node_destroy

    subroutine c_catalyst_conduit_node_set_path_external_float64_ptr(cnode, path, data, num_elements) 
        use iso_c_binding
        implicit none
        type(C_PTR), value, intent(IN) :: cnode
        character(kind=C_CHAR), intent(IN) :: path(*)
        real(8), intent (IN), dimension (*) :: data
        integer(C_SIZE_T), value, intent(in) :: num_elements
    end subroutine c_catalyst_conduit_node_set_path_external_float64_ptr

endmodule catalyst_conduit

module catalyst_api
    enum, bind(c)

        ! Use this with `integer(kind(catalyst_status))` to declare variables of this type.
        enumerator :: catalyst_status = 0

        enumerator :: catalyst_status_ok = 0
        enumerator :: catalyst_status_error_no_implementation = 1
        enumerator :: catalyst_status_error_already_loaded = 2
        enumerator :: catalyst_status_error_not_found = 3
        enumerator :: catalyst_status_error_not_catalyst = 4
        enumerator :: catalyst_status_error_incomplete = 5
        enumerator :: catalyst_status_error_unsupported_version = 6
        enumerator :: catalyst_status_error_conduit_mismatch = 7

    end enum

contains
    function c_catalyst_initialize(cnode) result(catalyst_code)
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none
        type(c_ptr), value, intent(in) :: cnode
        integer(c_int) :: catalyst_code
        catalyst_code = 0
    end function c_catalyst_initialize

    function c_catalyst_finalize(node) result(catalyst_code)
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none
        type(c_ptr), value, intent(in) :: node
        integer(c_int) :: catalyst_code
        catalyst_code = 0
    end function c_catalyst_finalize

    function c_catalyst_execute(node) result(catalyst_code)
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none
        type(c_ptr), value, intent(in) :: node
        integer(c_int) :: catalyst_code
        catalyst_code = 0
    end function c_catalyst_execute

    function c_catalyst_about(node) result(catalyst_code) 
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none
        type(c_ptr), value :: node
        integer(c_int) :: catalyst_code
        catalyst_code = 0
    end function c_catalyst_about

    function c_catalyst_results(node) result(catalyst_status) 
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none
        type(c_ptr), value :: node
        integer(c_int) :: catalyst_status
        catalyst_status = 0
    end function c_catalyst_results

endmodule catalyst_api
