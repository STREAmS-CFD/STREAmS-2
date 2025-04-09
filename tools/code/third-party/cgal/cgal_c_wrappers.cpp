#include "cgal_wrappers.h"

extern "C" {

    void polyhedron_from_file (Polytree **ptree, const char *fname){
        
  //    cout << "reading file: " << fname <<  std::endl;

        Polyhedron *P = new Polyhedron;

        std::ifstream input(fname, std::ios::in | std::ios::binary);

        //-----------------------------------------------------------
        // START - STL read
        //-----------------------------------------------------------
        //CGAL::read_STL( input,
        //    points,
        //    faces,
        //    true);
        //std::cout.precision(17);
        //std::cout << "OFF\n" << points.size() << " " << faces.size()  << " 0" << std::endl;
        //for(std::size_t i=0; i < points.size(); i++){
        //    std::cout << points[i][0] << " " << points[i][1] << " " << points[i][2]<< std::endl;
        //}
        //for(std::size_t i=0; i < faces.size(); i++){
        //    std::cout << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << std::endl;
        //}
        //-----------------------------------------------------------
        // END - STL read
        //-----------------------------------------------------------

        std::ifstream stream(fname);
        if(!stream) {
            std::cerr << "Cannot open file!" << std::endl;
        }
        stream >> *P;
        if(!stream) {
            std::cerr << "this .off file is not a polyhedron" << std::endl;
        }

        Tree *tree = new Tree(faces(*P).first, faces(*P).second, *P);
        tree->accelerate_distance_queries();

 //     std::cout << " n. facets: "    << P->size_of_facets()    << std::endl;
 //     std::cout << " n. halfedges: " << P->size_of_halfedges() << std::endl;
 //     std::cout << " n. vertices: "  << P->size_of_vertices()  << std::endl;

        *ptree = new Polytree;
        (*ptree)->poly = P;
        (*ptree)->tree = tree;

    }

    void polyhedron_closest (Polytree *ptree, REAL query_x, REAL query_y, REAL query_z,
            REAL *near_x, REAL *near_y, REAL *near_z){

        Point query_point(query_x,query_y,query_z);

        Point closest = ptree->tree->closest_point(query_point);

        *near_x = closest.x();
        *near_y = closest.y();
        *near_z = closest.z();
    }

    void polyhedron_closest_triangle (Polytree *ptree, REAL query_x, REAL query_y, REAL query_z,
            REAL *tri_1_x, REAL *tri_1_y, REAL *tri_1_z, 
            REAL *tri_2_x, REAL *tri_2_y, REAL *tri_2_z, 
            REAL *tri_3_x, REAL *tri_3_y, REAL *tri_3_z) {

        Point query_point(query_x,query_y,query_z);

        Point_and_primitive_id pp = ptree->tree->closest_point_and_primitive(query_point);
        Point closest_point = pp.first;
        Polyhedron::Face_handle f = pp.second; // closest primitive id
        *tri_1_x = f->halfedge()->vertex()->point().hx();
        *tri_1_y = f->halfedge()->vertex()->point().hy();
        *tri_1_z = f->halfedge()->vertex()->point().hz();
        *tri_2_x = f->halfedge()->next()->vertex()->point().hx();
        *tri_2_y = f->halfedge()->next()->vertex()->point().hy();
        *tri_2_z = f->halfedge()->next()->vertex()->point().hz();
        *tri_3_x = f->halfedge()->next()->next()->vertex()->point().hx();
        *tri_3_y = f->halfedge()->next()->next()->vertex()->point().hy();
        *tri_3_z = f->halfedge()->next()->next()->vertex()->point().hz();
        //std::cout << "closest point: " << closest_point << std::endl;
        //std::cout << "closest triangle: ( "
        //      << f->halfedge()->vertex()->point() << " , "
        //      << f->halfedge()->next()->vertex()->point() << " , "
        //      << f->halfedge()->next()->next()->vertex()->point()
        //      << " )" << std::endl;
    }

    bool polyhedron_inside(Polytree *ptree, REAL query_x, REAL query_y, REAL query_z) {

        //Point_inside inside_tester(*(ptree->tree));
        Point query_point = Point(query_x,query_y,query_z);

        //std::cout << "testing inside of point: " << query_x << " " << query_y << " " << query_z << std::endl;

        //RIMETTERECGAL::Side_of_triangle_mesh<Polyhedron, K> inside_tester(*(ptree->poly));
        CGAL::Side_of_triangle_mesh<Polyhedron, K> inside_tester(*(ptree->tree));

        // Determine the side and return true if inside!
        int ret=0;
        bool is_inside = (inside_tester(query_point) == CGAL::ON_BOUNDED_SIDE);
        if(is_inside)
            ret = 1;
        //std::cout << "returning inside: " << ret << std::endl;
        return ret;
    }

        void polyhedron_intersection(Polytree *ptree, REAL p1_x, REAL p1_y, REAL p1_z,
            REAL p2_x, REAL p2_y, REAL p2_z, REAL *inters_x, REAL *inters_y, REAL *inters_z, int* inters_type) {
        // constructs segment query
        Point a(p1_x, p1_y, p1_z);
        Point b(p2_x, p2_y, p2_z);
        Segment segment_query(a,b);
        // tests intersections with segment query
        if(ptree->tree->do_intersect(segment_query)) {
       //     std::cout << "intersection(s)" << std::endl;
        }
        else {
       //     std::cout << "no intersection" << std::endl;
            *inters_x = 0.0;
            *inters_y = 0.0;
            *inters_z = 0.0;
            *inters_type = 0;
        }

        // computes #intersections with segment query
        // std::cout << ptree->tree->number_of_intersected_primitives(segment_query) << " intersection(s)" << std::endl;
        // computes first encountered intersection with segment query (generally a point)
        Segment_intersection intersection = ptree->tree->any_intersection(segment_query);
        if(intersection) // gets intersection object
        {
            const Point* p = boost::get<Point>(&(intersection->first));
            if(p) {
       //         std::cout << "intersection object is a point " << *p << std::endl;
                *inters_x = p->hx();
                *inters_y = p->hy();
                *inters_z = p->hz();
                *inters_type = 1;
            } else {
                std::cout << "intersection is not a point " << *p << std::endl;
                *inters_x = 0.0;
                *inters_y = 0.0;
                *inters_z = 0.0;
                *inters_type = -1;
            }
        }
        // // computes all intersections with segment query (as pairs object - primitive_id)
        // std::list<Segment_intersection> intersections;
        // tree.all_intersections(segment_query, std::back_inserter(intersections));
        // // computes all intersected primitives with segment query as primitive ids
        // std::list<Primitive_id> primitives;
        // tree.all_intersected_primitives(segment_query, std::back_inserter(primitives));
    }

    void polyhedron_bbox(Polytree *ptree,
            REAL *bmin_x, REAL *bmin_y, REAL *bmin_z,
            REAL *bmax_x, REAL *bmax_y, REAL *bmax_z){
        Bbox_3 bbox = ptree->tree->bbox();
        *bmin_x = bbox.xmin();
        *bmin_y = bbox.ymin();
        *bmin_z = bbox.zmin();
        *bmax_x = bbox.xmax();
        *bmax_y = bbox.ymax();
        *bmax_z = bbox.zmax();
    }

    void polyhedron_finalize(Polytree **ptree){
        delete (*ptree)->tree; (*ptree)->tree = NULL;
        delete (*ptree)->poly; (*ptree)->poly = NULL;
        delete *ptree; *ptree = NULL;
    }

}
