#include <cstdlib>
#include <cstdio>

extern "C" {

    void polyhedron_from_file (void **ptree, const char *fname){
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    void polyhedron_closest (void *ptree, double query_x, double query_y, double query_z,
            double *near_x, double *near_y, double *near_z){
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    void polyhedron_closest_triangle (void *ptree, double query_x, double query_y, double query_z,
            double *tri_1_x, double *tri_1_y, double *tri_1_z, 
            double *tri_2_x, double *tri_2_y, double *tri_2_z, 
            double *tri_3_x, double *tri_3_y, double *tri_3_z) {
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    bool polyhedron_inside(void *ptree, double query_x, double query_y, double query_z) {
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
        return 0;
    }

    void polyhedron_intersection(void *ptree, double p1_x, double p1_y, double p1_z,
            double p2_x, double p2_y, double p2_z, double *inters_x, double *inters_y, double *inters_z, int* inters_type) {
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    void polyhedron_bbox(void *ptree,
            double *bmin_x, double *bmin_y, double *bmin_z,
            double *bmax_x, double *bmax_y, double *bmax_z){
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

    void polyhedron_finalize(void **ptree){
        printf("CGAL compilation not enabled! Exiting...\n");
        exit(0);
    }

}
