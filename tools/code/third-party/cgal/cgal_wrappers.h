#ifndef IMMERSED_H
#define IMMERSED_H

#ifdef SINGLE_PRECISION
#define REAL float
#else
#define REAL double
#endif

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

typedef CGAL::Simple_cartesian<REAL> K;

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

typedef struct {REAL x,y,z;} d3;

typedef struct {Polyhedron *poly; Tree *tree;} Polytree;

typedef boost::optional< Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;

extern "C" {
void polyhedron_from_file (Polytree  **ptree, const char *fname);

void polyhedron_closest (Polytree *ptree, REAL query_x, REAL query_y, REAL query_z,
                REAL *near_x, REAL *near_y, REAL *near_z);
void polyhedron_closest_triangle (Polytree *ptree, REAL query_x, REAL query_y, REAL query_z,
            REAL *tri_1_x, REAL *tri_1_y, REAL *tri_1_z, 
            REAL *tri_2_x, REAL *tri_2_y, REAL *tri_2_z, 
            REAL *tri_3_x, REAL *tri_3_y, REAL *tri_3_z);

bool polyhedron_inside(Polytree *ptree, REAL query_x, REAL query_y, REAL query_z);

void polyhedron_intersection(Polytree *ptree, REAL p1_x, REAL p1_y, REAL p1_z,
            REAL p2_x, REAL p2_y, REAL p2_z, REAL *inters_x, REAL *inters_y, REAL *inters_z, int* inters_type);

void polyhedron_bbox(Polytree *ptree, REAL *bmin_x, REAL *bmin_y, REAL *bmin_z, REAL *bmax_x, REAL *bmax_y, REAL *bmax_z);

void polyhedron_finalize(Polytree **ptree);
}

using std::cout;
using std::endl;

extern "C" int debuglevel;

#endif
