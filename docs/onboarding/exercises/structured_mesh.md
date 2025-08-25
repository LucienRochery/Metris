[Back to main.md](../main.md)

# Structured mesh

This function is implemented in [bunit/doc_structured_mesh.cxx](../../../bunit/doc_ball.hxx) and tested in [bunit/doc_structured_mesh.cxx](../../../bunit/doc_ball.cxx). 

Inputs:
- Bounds `xmin, xmax` and `ymin, ymax` 
- Grid counts `nx` and `ny`

Outputs: 
- Connectivity array of triangles (same as `fac2poi` in Metris)
- Connectivity array of edges
- Coordinates array (same as `coord`)
- Neighbour arrays for each (same as `fac2fac` and `edg2edg`)
- List of geometric nodes (corners)

Tips: 
- In the header file, write the function `doc_structured_mesh` which creates these arrays.
- In the unit test, call the function then feed the outputs to a `MetrisAPI` object, and then that to a `MetrisRunner` object, to get a `MeshBase` out of it. This is a valid mesh object that has, for instance, reconstructed its own neighbour arrays.
- Check this has the exepected number of points and elements. 
- Test the neighbour arrays by making sure neighbours are neighbours of each other and share a same facet. Functions in `aux_topo.hxx` such as `getedgfac` can be used for this. 



Notes:
- The function is implemented in a header file so it can later be used in other exercises, using the the `MetrisAPI` into `MetrisRunner`. 