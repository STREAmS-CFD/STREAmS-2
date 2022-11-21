module cgal_wrappers
    use iso_c_binding
    implicit none
    interface
        subroutine polyhedron_from_file(ptree, fname) bind(C, name="polyhedron_from_file")
            import
            type(c_ptr) :: ptree
            character(kind=c_char), intent(in) :: fname(*)
        end subroutine

        subroutine polyhedron_closest(ptree, query_x, query_y, query_z, near_x, near_y, near_z) bind(C,name="polyhedron_closest")
          import
          type(c_ptr),value :: ptree
          real(c_double), value :: query_x, query_y, query_z
          real(c_double) :: near_x, near_y, near_z
        end subroutine

        subroutine polyhedron_closest_triangle(ptree, query_x, query_y, query_z, &
            tri_1_x, tri_1_y, tri_1_z, &
            tri_2_x, tri_2_y, tri_2_z, &
            tri_3_x, tri_3_y, tri_3_z) &
            bind(C,name="polyhedron_closest_triangle")
          import
          type(c_ptr),value :: ptree
          real(c_double), value :: query_x, query_y, query_z
          real(c_double) :: tri_1_x, tri_1_y, tri_1_z
          real(c_double) :: tri_2_x, tri_2_y, tri_2_z
          real(c_double) :: tri_3_x, tri_3_y, tri_3_z
        end subroutine
        
        function polyhedron_inside(ptree, query_x, query_y, query_z) result(res) bind(C,name="polyhedron_inside")
          import
          logical(c_bool) :: res
          type(c_ptr),value :: ptree
          real(c_double), value :: query_x, query_y, query_z
        end function
        
        subroutine polyhedron_intersection(ptree, p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, &
            inters_x, inters_y, inters_z, inters_type) bind(C,name="polyhedron_intersection")
          import
          logical(c_bool) :: res
          type(c_ptr),value :: ptree
          real(c_double), value :: p1_x, p1_y, p1_z, p2_x, p2_y, p2_z
          real(c_double) :: inters_x, inters_y, inters_z
          integer(c_int) :: inters_type
        endsubroutine

        subroutine polyhedron_bbox(ptree, bmin_x, bmin_y, bmin_z, bmax_x, bmax_y, bmax_z) bind(C,name="polyhedron_bbox")
          import
          type(c_ptr),value :: ptree
          real(c_double) :: bmin_x, bmin_y, bmin_z, bmax_x, bmax_y, bmax_z
        end subroutine
    endinterface
contains
    subroutine cgal_polyhedron_read(ptree, fname)
        type(c_ptr),intent(out) :: ptree
        character(*),intent(in) :: fname

        call polyhedron_from_file(ptree, trim(adjustl(fname))//c_null_char)
    endsubroutine cgal_polyhedron_read

    function cgal_polyhedron_inside(ptree, query_x, query_y, query_z) result(is_inside)
        logical :: is_inside
        type(c_ptr),value :: ptree
        real(c_double) :: query_x, query_y, query_z

        is_inside = polyhedron_inside(ptree, query_x, query_y, query_z)
        !print*,'is_inside wrapper: ',is_inside
    endfunction cgal_polyhedron_inside

        subroutine cgal_polyhedron_intersection(ptree, p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, inters_x, &
                   inters_y, inters_z, inters_type)
        type(c_ptr),value :: ptree
        real(c_double) :: p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, inters_x, inters_y, inters_z
        integer(c_int) :: inters_type

        call polyhedron_intersection(ptree, p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, inters_x, inters_y, &
                inters_z, inters_type)
        !print*,'is_inside wrapper: ',is_inside
    endsubroutine cgal_polyhedron_intersection
endmodule cgal_wrappers
