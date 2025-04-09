#include <cstdlib>
#include <cstdio>

#ifdef SINGLE_PRECISION
#define REAL float
#else
#define REAL double
#endif

extern "C" {

    void polyhedron_from_file (void **ptree, const char *fname){
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    void polyhedron_closest (void *ptree, REAL query_x, REAL query_y, REAL query_z,
            REAL *near_x, REAL *near_y, REAL *near_z){
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    void polyhedron_closest_triangle (void *ptree, REAL query_x, REAL query_y, REAL query_z,
            REAL *tri_1_x, REAL *tri_1_y, REAL *tri_1_z, 
            REAL *tri_2_x, REAL *tri_2_y, REAL *tri_2_z, 
            REAL *tri_3_x, REAL *tri_3_y, REAL *tri_3_z) {
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    bool polyhedron_inside(void *ptree, REAL query_x, REAL query_y, REAL query_z) {
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
        return 0;
    }

    void polyhedron_intersection(void *ptree, REAL p1_x, REAL p1_y, REAL p1_z,
            REAL p2_x, REAL p2_y, REAL p2_z, REAL *inters_x, REAL *inters_y, REAL *inters_z, int* inters_type) {
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    void polyhedron_bbox(void *ptree,
            REAL *bmin_x, REAL *bmin_y, REAL *bmin_z,
            REAL *bmax_x, REAL *bmax_y, REAL *bmax_z){
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    void polyhedron_finalize(void **ptree){
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

}
