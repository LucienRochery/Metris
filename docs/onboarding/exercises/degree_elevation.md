[Back to main.md](../main.md)

# Degree elevation


This function is implemented in [bunit/doc_degree_elevation.cxx](../../../bunit/doc_degree_elevation.cxx).

Inputs:
- Element connectivity arrays `edg2poi`, `fac2poi`, `tet2poi`
- Coordinates array `coord`
- A target degree `tardeg`
- Possibly tag arrays (see [ball exercise](ball.md))

Outputs:
- Arrays now represent a degree `tardeg` mesh. This must be conforming, i.e. two elements sharing a facet must also share its control points. 
- The output mesh is straight, i.e. control points (or Lagrange nodes, equivalently since straight) are uniformly spaced. 


Notes:
- A first implementation can be done that only considers the highest dimension elements (i.e. triangles in 2D and tetrahedra in 3D). 
- `tardeg = 2` only can be targetted at first. For higher degrees, it is necessary to compile Metris using `cmake -DMETRIS_MAX_DEG=<n>`, but this is possibly currently disabled on the GitHub version of Metris. 

Tips:
- Coordinates array can be expanded using `inc_n()`. This takes care of allocations if necessary. 
- See the [high-order element evaluation exercise](HO_evaluation.md) for HO ordering tips. 
- Connectivity arrays can be resized using `allocate(nelem, nnode)`, this modifies the stride (second dimension) to match `nnode`.
- The main difficulty of this exercise is dealing with node indices and making sure facets aren't mangled when copied from one element to another. A helper function can be written for this. 

Testing:
- A simple way to check the routine is working properly is writing out the mesh using the `writeMesh` function that takes loose arrays (`io_libmeshb.hxx`)

