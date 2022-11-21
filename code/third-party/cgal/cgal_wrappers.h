#ifndef IMMERSED_H
#define IMMERSED_H

// Author(s) : Pierre Alliez
#include <fstream>
#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Side_of_triangle_mesh.h>

typedef CGAL::Simple_cartesian<double> K;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

typedef CGAL::Side_of_triangle_mesh<Polyhedron, K> Point_inside;
typedef CGAL::Bbox_3 Bbox_3;

typedef struct {double x,y,z;} d3;

typedef struct {Polyhedron *poly; Tree *tree;} Polytree;

typedef boost::optional< Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;

extern "C" {
void polyhedron_from_file (Polytree  **ptree, const char *fname);

void polyhedron_closest (Polytree *ptree, double query_x, double query_y, double query_z,
                double *near_x, double *near_y, double *near_z);
void polyhedron_closest_triangle (Polytree *ptree, double query_x, double query_y, double query_z,
            double *tri_1_x, double *tri_1_y, double *tri_1_z, 
            double *tri_2_x, double *tri_2_y, double *tri_2_z, 
            double *tri_3_x, double *tri_3_y, double *tri_3_z);

bool polyhedron_inside(Polytree *ptree, double query_x, double query_y, double query_z);

void polyhedron_intersection(Polytree *ptree, double p1_x, double p1_y, double p1_z,
            double p2_x, double p2_y, double p2_z, double *inters_x, double *inters_y, double *inters_z, int* inters_type);

void polyhedron_bbox(Polytree *ptree, double *bmin_x, double *bmin_y, double *bmin_z, double *bmax_x, double *bmax_y, double *bmax_z);

void polyhedron_finalize(Polytree **ptree);
}

using std::cout;
using std::endl;

extern "C" int debuglevel;

#endif
