/* kpoint.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include "mathfunc.h"
#include "symmetry.h"
#include "kpoint.h"

#include "debug.h"

/* #define GRID_ORDER_XYZ */
/* The addressing order of mesh grid is defined as running left */
/* element first. But when GRID_ORDER_XYZ is defined, it is changed to right */ 
/* element first. */

static int search_space[][3] = {{0, 0, 0},
				{0, 0, 1},
				{0, 1, -1},
				{0, 1, 0},
				{0, 1, 1},
				{1, -1, -1},
				{1, -1, 0},
				{1, -1, 1},
				{1, 0, -1},
				{1, 0, 0},
				{1, 0, 1},
				{1, 1, -1},
				{1, 1, 0},
				{1, 1, 1},
				{-1, -1, -1},
				{-1, -1, 0},
				{-1, -1, 1},
				{-1, 0, -1},
				{-1, 0, 0},
				{-1, 0, 1},
				{-1, 1, -1},
				{-1, 1, 0},
				{-1, 1, 1},
				{0, -1, -1},
				{0, -1, 0},
				{0, -1, 1},
				{0, 0, -1}};

static PointSymmetry get_point_group_reciprocal(const MatINT * rotations,
						const int is_time_reversal);
static PointSymmetry
get_point_group_reciprocal_with_q(SPGCONST PointSymmetry * pointgroup,
				  const double symprec,
				  const int num_q,
				  SPGCONST double qpoints[][3]);
static int get_ir_kpoints(int map[],
			  SPGCONST double kpoints[][3],
			  const int num_kpoint,
			  SPGCONST PointSymmetry * point_symmetry,
			  const double symprec);
static int get_ir_reciprocal_mesh(int grid_address[][3],
				  int map[],
				  const int mesh[3],
				  const int is_shift[3],
				  SPGCONST PointSymmetry * point_symmetry);
static int
get_ir_reciprocal_mesh_openmp(int grid_address[][3],
			      int map[],
			      const int mesh[3],
			      const int is_shift[3],
			      SPGCONST PointSymmetry * point_symmetry);
static int relocate_BZ_grid_address(int bz_grid_address[][3],
				    int bz_map[],
				    int grid_address[][3],
				    const int mesh[3],
				    SPGCONST double rec_lattice[3][3],
				    const int is_shift[3]);
static double get_tolerance_for_BZ_reduction(SPGCONST double rec_lattice[3][3]);
static int get_ir_triplets_at_q(int weights[],
				int grid_address[][3],
				int third_q[],
				const int grid_point,
				const int mesh[3],
				SPGCONST PointSymmetry * pointgroup);
static int get_BZ_triplets_at_q(int triplets[][3],
				const int grid_point,
				SPGCONST int bz_grid_address[][3],
				const int bz_map[],
				const int weights[],
				const int mesh[3]);
static void get_third_q_of_triplets_at_q(int address[3][3],
					 const int bz_map[],
					 const int mesh[3],
					 const int bzmesh[3],
					 const int bzmesh_double[3]);
static void grid_point_to_grid_double(int grid_double[3],
				      const int address,
				      const int mesh[3],
				      const int is_shift[3]);
static void get_grid_address(int grid_address[3],
			     const int grid[3],
			     const int mesh[3]);
static void get_vector_modulo(int v[3], const int m[3]);
static int get_grid_point(const int grid[3],
			  const int mesh[3]);

int kpt_get_irreducible_kpoints(int map[],
				SPGCONST double kpoints[][3],
				const int num_kpoint,
				const Symmetry * symmetry,
				const int is_time_reversal,
				const double symprec)
{
  int i;
  PointSymmetry point_symmetry;
  MatINT *rotations;
  
  rotations = mat_alloc_MatINT(symmetry->size);
  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(rotations->mat[i], symmetry->rot[i]);
  }

  point_symmetry = get_point_group_reciprocal(rotations,
					      is_time_reversal);
  mat_free_MatINT(rotations);

  return get_ir_kpoints(map, kpoints, num_kpoint, &point_symmetry, symprec);
}

/* grid_address (e.g. 4x4x4 mesh)                             */
/*    [[ 0  0  0]                                             */
/*     [ 1  0  0]                                             */
/*     [ 2  0  0]                                             */
/*     [-1  0  0]                                             */
/*     [ 0  1  0]                                             */
/*     [ 1  1  0]                                             */
/*     [ 2  1  0]                                             */
/*     [-1  1  0]                                             */
/*     ....      ]                                            */
/*                                                            */
/* Each value of 'map' correspnds to the index of grid_point. */
int kpt_get_irreducible_reciprocal_mesh(int grid_address[][3],
					int map[],
					const int mesh[3],
					const int is_shift[3],
					const int is_time_reversal,
					const Symmetry * symmetry)
{
  int i;
  PointSymmetry point_symmetry;
  MatINT *rotations;
  
  rotations = mat_alloc_MatINT(symmetry->size);
  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(rotations->mat[i], symmetry->rot[i]);
  }

  point_symmetry = get_point_group_reciprocal(rotations,
					      is_time_reversal);
  mat_free_MatINT(rotations);

#ifdef _OPENMP
  return get_ir_reciprocal_mesh_openmp(grid_address,
				       map,
				       mesh,
				       is_shift,
				       &point_symmetry);
#else
  return get_ir_reciprocal_mesh(grid_address,
				map,
				mesh,
				is_shift,
				&point_symmetry);
#endif
  
}

int kpt_get_stabilized_reciprocal_mesh(int grid_address[][3],
				       int map[],
				       const int mesh[3],
				       const int is_shift[3],
				       const int is_time_reversal,
				       const MatINT * rotations,
				       const int num_q,
				       SPGCONST double qpoints[][3])
{
  PointSymmetry pointgroup, pointgroup_q;
  double tolerance;
  
  pointgroup = get_point_group_reciprocal(rotations,
					  is_time_reversal);

  tolerance = 0.1 / (mesh[0] + mesh[1] + mesh[2]);
  pointgroup_q = get_point_group_reciprocal_with_q(&pointgroup,
						   tolerance,
						   num_q,
						   qpoints);

#ifdef _OPENMP
  return get_ir_reciprocal_mesh_openmp(grid_address,
				       map,
				       mesh,
				       is_shift,
				       &pointgroup_q);
#else
  return get_ir_reciprocal_mesh(grid_address,
				map,
				mesh,
				is_shift,
				&pointgroup_q);
#endif

}

int kpt_relocate_BZ_grid_address(int bz_grid_address[][3],
				 int bz_map[],
				 int grid_address[][3],
				 const int mesh[3],
				 SPGCONST double rec_lattice[3][3],
				 const int is_shift[3])
{
  return relocate_BZ_grid_address(bz_grid_address,
				  bz_map,
				  grid_address,
				  mesh,
				  rec_lattice,
				  is_shift);
}

int kpt_get_ir_triplets_at_q(int weights[],
			     int grid_address[][3],
			     int third_q[],
			     const int grid_point,
			     const int mesh[3],
			     const int is_time_reversal,
			     const MatINT * rotations)
{
  PointSymmetry pointgroup;

  pointgroup = get_point_group_reciprocal(rotations,
					  is_time_reversal);
  return get_ir_triplets_at_q(weights,
			      grid_address,
			      third_q,
			      grid_point,
			      mesh,
			      &pointgroup);
}

int kpt_get_BZ_triplets_at_q(int triplets[][3],
			     const int grid_point,
			     SPGCONST int bz_grid_address[][3],
			     const int bz_map[],
			     const int weights[],
			     const int mesh[3])
{
  return get_BZ_triplets_at_q(triplets,
			      grid_point,
			      bz_grid_address,
			      bz_map,
			      weights,
			      mesh);
}


      

static PointSymmetry get_point_group_reciprocal(const MatINT * rotations,
						const int is_time_reversal)
{
  int i, j, num_pt = 0;
  MatINT *rot_reciprocal;
  PointSymmetry point_symmetry;
  SPGCONST int inversion[3][3] = {
    {-1, 0, 0 },
    { 0,-1, 0 },
    { 0, 0,-1 }
  };
  
  if (is_time_reversal) {
    rot_reciprocal = mat_alloc_MatINT(rotations->size * 2);
  } else {
    rot_reciprocal = mat_alloc_MatINT(rotations->size);
  }

  for (i = 0; i < rotations->size; i++) {
    mat_transpose_matrix_i3(rot_reciprocal->mat[i], rotations->mat[i]);
    
    if (is_time_reversal) {
      mat_multiply_matrix_i3(rot_reciprocal->mat[rotations->size+i],
			     inversion,
			     rot_reciprocal->mat[i]);
    }
  }


  for (i = 0; i < rot_reciprocal->size; i++) {
    for (j = 0; j < num_pt; j++) {
      if (mat_check_identity_matrix_i3(point_symmetry.rot[j],
				       rot_reciprocal->mat[i])) {
	goto escape;
      }
    }
    
    mat_copy_matrix_i3(point_symmetry.rot[num_pt],
		       rot_reciprocal->mat[i]);
    num_pt++;
  escape:
    ;
  }

  point_symmetry.size = num_pt;

  mat_free_MatINT(rot_reciprocal);

  return point_symmetry;
}

static PointSymmetry
get_point_group_reciprocal_with_q(SPGCONST PointSymmetry * pointgroup,
				  const double symprec,
				  const int num_q,
				  SPGCONST double qpoints[][3])
{
  int i, j, k, l, is_all_ok=0, num_ptq = 0;
  double q_rot[3], diff[3];
  PointSymmetry pointgroup_q;

  for (i = 0; i < pointgroup->size; i++) {
    for (j = 0; j < num_q; j++) {
      is_all_ok = 0;
      mat_multiply_matrix_vector_id3(q_rot,
				     pointgroup->rot[i],
				     qpoints[j]);

      for (k = 0; k < num_q; k++) {
	for (l = 0; l < 3; l++) {
	  diff[l] = q_rot[l] - qpoints[k][l];
	  diff[l] -= mat_Nint(diff[l]);
	}
	
	if (mat_Dabs(diff[0]) < symprec && 
	    mat_Dabs(diff[1]) < symprec &&
	    mat_Dabs(diff[2]) < symprec) {
	  is_all_ok = 1;
	  break;
	}
      }
      
      if (! is_all_ok) {
	break;
      }
    }

    if (is_all_ok) {
      mat_copy_matrix_i3(pointgroup_q.rot[num_ptq], pointgroup->rot[i]);
      num_ptq++;
    }
  }
  pointgroup_q.size = num_ptq;

  return pointgroup_q;
}


static int get_ir_kpoints(int map[],
			  SPGCONST double kpoints[][3],
			  const int num_kpoint,
			  SPGCONST PointSymmetry * point_symmetry,
			  const double symprec)
{
  int i, j, k, l, num_ir_kpoint = 0, is_found;
  int *ir_map;
  double kpt_rot[3], diff[3];

  ir_map = (int*)malloc(num_kpoint*sizeof(int));

  for (i = 0; i < num_kpoint; i++) {

    map[i] = i;

    is_found = 1;

    for (j = 0; j < point_symmetry->size; j++) {
      mat_multiply_matrix_vector_id3(kpt_rot, point_symmetry->rot[j], kpoints[i]);

      for (k = 0; k < 3; k++) {
	diff[k] = kpt_rot[k] - kpoints[i][k];
	diff[k] = diff[k] - mat_Nint(diff[k]);
      }

      if (mat_Dabs(diff[0]) < symprec && 
	  mat_Dabs(diff[1]) < symprec && 
	  mat_Dabs(diff[2]) < symprec) {
	continue;
      }
      
      for (k = 0; k < num_ir_kpoint; k++) {
	mat_multiply_matrix_vector_id3(kpt_rot, point_symmetry->rot[j], kpoints[i]);

	for (l = 0; l < 3; l++) {
	  diff[l] = kpt_rot[l] - kpoints[ir_map[k]][l];
	  diff[l] = diff[l] - mat_Nint(diff[l]);
	}

	if (mat_Dabs(diff[0]) < symprec && 
	    mat_Dabs(diff[1]) < symprec && 
	    mat_Dabs(diff[2]) < symprec) {
	  is_found = 0;
	  map[i] = ir_map[k];
	  break;
	}
      }

      if (! is_found)
	break;
    }

    if (is_found) {
      ir_map[num_ir_kpoint] = i;
      num_ir_kpoint++;
    }
  }

  free(ir_map);
  ir_map = NULL;

  return num_ir_kpoint;
}

static int get_ir_reciprocal_mesh(int grid_address[][3],
				  int map[],
				  const int mesh[3],
				  const int is_shift[3],
				  SPGCONST PointSymmetry * point_symmetry)
{
  /* In the following loop, mesh is doubled. */
  /* Even and odd mesh numbers correspond to */
  /* is_shift[i] = 0 and 1, respectively. */
  /* is_shift = [0,0,0] gives Gamma center mesh. */
  /* grid: reducible grid points */
  /* map: the mapping from each point to ir-point. */

  int i, j, k, l, grid_point, grid_point_rot, num_ir = 0;
  int grid_double[3], grid_rot[3], mesh_double[3];

  for (i = 0; i < 3; i++) {
    mesh_double[i] = mesh[i] * 2;
  }

  /* "-1" means the element is not touched yet. */
  for (i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++) {
    map[i] = -1;
  }

#ifndef GRID_ORDER_XYZ
  for (i = 0; i < mesh[2]; i++) {
    for (j = 0; j < mesh[1]; j++) {
      for (k = 0; k < mesh[0]; k++) {
	grid_double[0] = k * 2 + is_shift[0];
	grid_double[1] = j * 2 + is_shift[1];
	grid_double[2] = i * 2 + is_shift[2];
#else
  for (i = 0; i < mesh[0]; i++) {
    for (j = 0; j < mesh[1]; j++) {
      for (k = 0; k < mesh[2]; k++) {
  	grid_double[0] = i * 2 + is_shift[0];
  	grid_double[1] = j * 2 + is_shift[1];
  	grid_double[2] = k * 2 + is_shift[2];
#endif	

	grid_point = get_grid_point(grid_double, mesh);
	get_grid_address(grid_address[grid_point], grid_double, mesh);

	for (l = 0; l < point_symmetry->size; l++) {
	  mat_multiply_matrix_vector_i3(grid_rot,
					point_symmetry->rot[l],	grid_double);
	  get_vector_modulo(grid_rot, mesh_double);
	  grid_point_rot = get_grid_point(grid_rot, mesh);

	  if (grid_point_rot > -1) { /* Invalid if even --> odd or odd --> even */
	    if (map[grid_point_rot] > -1) {
	      map[grid_point] = map[grid_point_rot];
	      break;
	    }
	  }
	}
	
	if (map[grid_point] == -1) {
	  map[grid_point] = grid_point;
	  num_ir++;
	}
      }
    }
  }

  return num_ir;
}

static int
get_ir_reciprocal_mesh_openmp(int grid_address[][3],
			      int map[],
			      const int mesh[3],
			      const int is_shift[3],
			      SPGCONST PointSymmetry * point_symmetry)
{
  int i, j, k, l, grid_point, grid_point_rot, num_ir;
  int grid_double[3], grid_rot[3], mesh_double[3];

  for (i = 0; i < 3; i++) {
    mesh_double[i] = mesh[i] * 2;
  }

#ifndef GRID_ORDER_XYZ
#pragma omp parallel for private(j, k, l, grid_point, grid_point_rot, grid_double, grid_rot)
  for (i = 0; i < mesh[2]; i++) {
    for (j = 0; j < mesh[1]; j++) {
      for (k = 0; k < mesh[0]; k++) {
	grid_double[0] = k * 2 + is_shift[0];
	grid_double[1] = j * 2 + is_shift[1];
	grid_double[2] = i * 2 + is_shift[2];
#else
#pragma omp parallel for private(j, k, l, grid_point, grid_point_rot, grid_double, grid_rot)
  for (i = 0; i < mesh[0]; i++) {
    for (j = 0; j < mesh[1]; j++) {
      for (k = 0; k < mesh[2]; k++) {
  	grid_double[0] = i * 2 + is_shift[0];
  	grid_double[1] = j * 2 + is_shift[1];
  	grid_double[2] = k * 2 + is_shift[2];
#endif	

	grid_point = get_grid_point(grid_double, mesh);
	map[grid_point] = grid_point;
	get_grid_address(grid_address[grid_point], grid_double, mesh);

	for (l = 0; l < point_symmetry->size; l++) {
	  mat_multiply_matrix_vector_i3(grid_rot,
					point_symmetry->rot[l],	grid_double);
	  get_vector_modulo(grid_rot, mesh_double);
	  grid_point_rot = get_grid_point(grid_rot, mesh);

	  if (grid_point_rot > -1) { /* Invalid if even --> odd or odd --> even */
	    if (grid_point_rot < map[grid_point]) {
	      map[grid_point] = grid_point_rot;
	    }
	  }
	}
      }
    }
  }

  num_ir = 0;

#pragma omp parallel for reduction(+:num_ir)
  for (i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++) {
    if (map[i] == i) {
      num_ir++;
    }
  }
  
  return num_ir;
}

/* Relocate grid addresses to first Brillouin zone */
/* bz_grid_address[prod(mesh + 1)][3] */
/* bz_map[prod(mesh * 2 - 1)] */
static int relocate_BZ_grid_address(int bz_grid_address[][3],
				    int bz_map[],
				    int grid_address[][3],
				    const int mesh[3],
				    SPGCONST double rec_lattice[3][3],
				    const int is_shift[3])
{
  double tolerance, min_distance;
  double vector[3], distance[27];
  int bzmesh[3], bzmesh_double[3], address_double[3];
  int i, j, k, min_index, boundary_gp, total_num_gp, bzgp, gp;

  tolerance = get_tolerance_for_BZ_reduction(rec_lattice);
  for (i = 0; i < 3; i++) {
    bzmesh[i] = mesh[i] * 2 - 1;
    bzmesh_double[i] = bzmesh[i] * 2;
  }
  for (i = 0; i < bzmesh[0] * bzmesh[1] * bzmesh[2]; i++) {
    bz_map[i] = -1;
  }
  
  boundary_gp = 0;
  total_num_gp = mesh[0] * mesh[1] * mesh[2];
  for (i = 0; i < total_num_gp; i++) {
    for (j = 0; j < 27; j++) {
      for (k = 0; k < 3; k++) {
	address_double[k] =
	  (grid_address[i][k] + search_space[j][k] * mesh[k]) * 2 + is_shift[k];
      }
      mat_multiply_matrix_vector_di3(vector, rec_lattice, address_double);
      distance[j] = mat_norm_squared_d3(vector);
    }
    min_distance = distance[0];
    min_index = 0;
    for (j = 1; j < 27; j++) {
      if (distance[j] + tolerance < min_distance) {
	min_distance = distance[j];
	min_index = j;
      }
    }

    for (j = 0; j < 27; j++) {
      if (distance[j] < min_distance + tolerance) {
	if (j == min_index) {
	  gp = i;
	} else {
	  gp = boundary_gp + total_num_gp;
	}
	for (k = 0; k < 3; k++) {
	  bz_grid_address[gp][k] = 
	    grid_address[i][k] + search_space[j][k] * mesh[k];
	  address_double[k] = bz_grid_address[gp][k] * 2 + is_shift[k];
	  if (address_double[k] < 0) {
	    address_double[k] += bzmesh_double[k];
	  }
	}
	bzgp = get_grid_point(address_double, bzmesh);
	bz_map[bzgp] = gp;
	if (j != min_index) {
	  boundary_gp++;
	}
      }
    }
  }

  return boundary_gp + total_num_gp;
}

static double get_tolerance_for_BZ_reduction(SPGCONST double rec_lattice[3][3])
{
  int i, j;
  double tolerance;
  double length[3];
  
  for (i = 0; i < 3; i++) {
    length[i] = 0;
    for (j = 0; j < 3; j++) {
      length[i] += rec_lattice[j][i] * rec_lattice[j][i];
    }
  }
  tolerance = length[0];
  for (i = 1; i < 3; i++) {
    if (tolerance > length[i]) {
      tolerance = length[i];
    }
  }
  tolerance /= 100;
  return tolerance;
}
 
static int get_ir_triplets_at_q(int weights[],
				int grid_address[][3],
				int third_q[],
				const int grid_point,
				const int mesh[3],
				SPGCONST PointSymmetry * pointgroup)
{
  int i, j, num_grid, q_2, num_ir_q, num_ir_triplets, ir_grid_point;
  int mesh_double[3], is_shift[3];
  int grid_double0[3], grid_double1[3], grid_double2[3];
  int *map_q, *ir_grid_points, *weight_q;
  double tolerance;
  double stabilizer_q[1][3];
  PointSymmetry pointgroup_q;

  tolerance = 0.1 / (mesh[0] + mesh[1] + mesh[2]);

  num_grid = mesh[0] * mesh[1] * mesh[2];

  for (i = 0; i < 3; i++) {
    /* Only consider the gamma-point */
    is_shift[i] = 0;
    mesh_double[i] = mesh[i] * 2;
  }

  /* Search irreducible q-points (map_q) with a stabilizer */
  grid_point_to_grid_double(grid_double0, grid_point, mesh, is_shift); /* q */
  for (i = 0; i < 3; i++) {
    stabilizer_q[0][i] =
      (double)grid_double0[i] / mesh_double[i] - (grid_double0[i] > mesh[i]);
  }

  pointgroup_q = get_point_group_reciprocal_with_q(pointgroup,
						   tolerance,
						   1,
						   stabilizer_q);
  map_q = (int*) malloc(sizeof(int) * num_grid);

#ifdef _OPENMP
  num_ir_q = get_ir_reciprocal_mesh_openmp(grid_address,
					   map_q,
					   mesh,
					   is_shift,
					   &pointgroup_q);
#else
  num_ir_q = get_ir_reciprocal_mesh(grid_address,
				    map_q,
				    mesh,
				    is_shift,
				    &pointgroup_q);
#endif

  ir_grid_points = (int*) malloc(sizeof(int) * num_ir_q);
  weight_q = (int*) malloc(sizeof(int) * num_grid);
  num_ir_q = 0;
  for (i = 0; i < num_grid; i++) {
    if (map_q[i] == i) {
      ir_grid_points[num_ir_q] = i;
      num_ir_q++;
    }
    weight_q[i] = 0;
    third_q[i] = -1;
    weights[i] = 0;
  }

  for (i = 0; i < num_grid; i++) {
    weight_q[map_q[i]]++;
  }

#pragma omp parallel for private(j, grid_double1, grid_double2)
  for (i = 0; i < num_ir_q; i++) {
    grid_point_to_grid_double(grid_double1, ir_grid_points[i], mesh, is_shift); /* q' */
    for (j = 0; j < 3; j++) { /* q'' */
      grid_double2[j] = - grid_double0[j] - grid_double1[j];
    }
    get_vector_modulo(grid_double2, mesh_double);
    third_q[ir_grid_points[i]] = get_grid_point(grid_double2, mesh);
  }

  num_ir_triplets = 0;
  for (i = 0; i < num_ir_q; i++) {
    ir_grid_point = ir_grid_points[i];
    q_2 = third_q[ir_grid_point];
    if (weights[map_q[q_2]]) {
      weights[map_q[q_2]] += weight_q[ir_grid_point];
    } else {
      weights[ir_grid_point] = weight_q[ir_grid_point];
      num_ir_triplets++;
    }
  }

  free(map_q);
  map_q = NULL;
  free(weight_q);
  weight_q = NULL;
  free(ir_grid_points);
  ir_grid_points = NULL;

  return num_ir_triplets;
}

static int get_BZ_triplets_at_q(int triplets[][3],
				const int grid_point,
				SPGCONST int bz_grid_address[][3],
				const int bz_map[],
				const int weights[],
				const int mesh[3])
{
  int i, j, k, num_ir;
  int address[3][3], address_double[3], bzmesh[3], bzmesh_double[3];
  int *ir_grid_points;

  for (i = 0; i < 3; i++) {
    bzmesh[i] = mesh[i] * 2 - 1;
    bzmesh_double[i] = bzmesh[i] * 2;
  }

  num_ir = 0;
  ir_grid_points = (int*) malloc(sizeof(int) * mesh[0] * mesh[1] * mesh[2]);
  for (i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++) {
    if (weights[i] > 0) {
      ir_grid_points[num_ir] = i;
      num_ir++;
    }
  }
 
#pragma omp parallel for private(j, k, address, address_double)
  for (i = 0; i < num_ir; i++) {
    for (j = 0; j < 3; j++) {
      address[0][j] = bz_grid_address[grid_point][j];
      address[1][j] = bz_grid_address[ir_grid_points[i]][j];
      address[2][j] = - address[0][j] - address[1][j];
    }
    get_third_q_of_triplets_at_q(address,
				 bz_map,
				 mesh,
				 bzmesh,
				 bzmesh_double);
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	address_double[k] = address[j][k] * 2;
	if (address_double[k] < 0) {
	  address_double[k] += bzmesh_double[k];
	}
      }
      triplets[i][j] = bz_map[get_grid_point(address_double, bzmesh)];
    }
  }

  free(ir_grid_points);
  
  return num_ir;
}

static void get_third_q_of_triplets_at_q(int address[3][3],
					 const int bz_map[],
					 const int mesh[3],
					 const int bzmesh[3],
					 const int bzmesh_double[3])
{
  int i, j, smallest_g, smallest_index, sum_g, delta_g[3];
  int bzgp[27], address_double[3];

  get_vector_modulo(address[2], mesh);
  for (i = 0; i < 3; i++) {
    delta_g[i] = 0;
    for (j = 0; j < 3; j++) {
      delta_g[i] += address[j][i];
    }
    delta_g[i] /= mesh[i];
  }
  
  for (i = 0; i < 27; i++) {
    for (j = 0; j < 3; j++) {
      address_double[j] = (address[2][j] + search_space[i][j] * mesh[j]) * 2;
    }
    if (abs(address_double[0] > bzmesh[0]) ||
	abs(address_double[1] > bzmesh[1]) ||
	abs(address_double[2] > bzmesh[2]) ||
	abs(address_double[0] < -bzmesh[0]) ||
	abs(address_double[1] < -bzmesh[1]) ||
	abs(address_double[2] < -bzmesh[2])) { /* outside extended zone */
      bzgp[i] = -1;
      continue;
    }
    for (j = 0; j < 3; j++) {
      if (address_double[j] < 0) {
	address_double[j] += bzmesh_double[j];
      }
    }
    bzgp[i] = bz_map[get_grid_point(address_double, bzmesh)];
  }

  for (i = 0; i < 27; i++) {
    if (bzgp[i] != -1) {
      goto escape;
    }
  }
  printf("******* Warning *******\n");
  printf(" No third-q was found.\n");
  printf("******* Warning *******\n");

 escape:

  smallest_g = 4;
  smallest_index = 0;
  for (i = 0; i < 27; i++) {
    if (bzgp[i] > -1) { /* q'' is in BZ */
      sum_g = (abs(delta_g[0] + search_space[i][0]) +
	       abs(delta_g[1] + search_space[i][1]) +
	       abs(delta_g[2] + search_space[i][2]));
      if (sum_g < smallest_g) {
	smallest_index = i;
	smallest_g = sum_g;
      }
    }
  }
  
  for (i = 0; i < 3; i++) {
    address[2][i] += search_space[smallest_index][i] * mesh[i];
  }
}

static int get_grid_point(const int grid_double[3],
			  const int mesh[3])
{
  int i, grid[3];

  for (i = 0; i < 3; i++) {
    if (grid_double[i] % 2 == 0) {
      grid[i] = grid_double[i] / 2;
    } else {
      grid[i] = (grid_double[i] - 1) / 2;
    }
  }

#ifndef GRID_ORDER_XYZ
  return grid[2] * mesh[0] * mesh[1] + grid[1] * mesh[0] + grid[0];
#else
  return grid[0] * mesh[1] * mesh[2] + grid[1] * mesh[2] + grid[2];
#endif  
}

static void grid_point_to_grid_double(int grid_double[3],
				      const int grid_point,
				      const int mesh[3],
				      const int is_shift[3])
{
  int i;
  int grid[3];

#ifndef GRID_ORDER_XYZ
  grid[2] = grid_point / (mesh[0] * mesh[1]);
  grid[1] = (grid_point - grid[2] * mesh[0] * mesh[1]) / mesh[0];
  grid[0] = grid_point % mesh[0];
#else
  grid[0] = grid_point / (mesh[1] * mesh[2]);
  grid[1] = (grid_point - grid[0] * mesh[1] * mesh[2]) / mesh[2];
  grid[2] = grid_point % mesh[2];
#endif

  for (i = 0; i < 3; i++) {
    grid_double[i] = grid[i] * 2 + is_shift[i];
  }
}

static void get_grid_address(int address[3],
			     const int grid_double[3],
			     const int mesh[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    if (grid_double[i] % 2 == 0) {
      address[i] = grid_double[i] / 2;
    } else {
      address[i] = (grid_double[i] - 1) / 2;
    }

#ifndef GRID_BOUNDARY_AS_NEGATIVE
    address[i] = address[i] - mesh[i] * (address[i] > mesh[i] / 2);
#else
    address[i] = address[i] - mesh[i] * (address[i] >= mesh[i] / 2);
#endif
  }  
}

static void get_vector_modulo(int v[3], const int m[3])
{
  int i;

  for (i = 0; i < 3; i++) {
    v[i] = v[i] % m[i];

    if (v[i] < 0)
      v[i] += m[i];
  }
}
