[Back to main.md](../main.md)

# Structured mesh

A function that produces the mesh data for a simple structured mesh can be useful in some validations.

This function is implemented in [bunit/doc_structured_mesh.cxx](../../../bunit/doc_ball.hxx) and tested in [bunit/doc_structured_mesh.cxx](../../../bunit/doc_ball.cxx).

Inputs:
- Bounds `xmin, xmax` and `ymin, ymax`
- Grid counts `nx` and `ny`

Outputs:
- Connectivity array of triangles (same as `fac2poi` in Metris)
- Connectivity array of edges
- Coordinates array (same as `coord`)
- Neighbor arrays for each (same as `fac2fac` and `edg2edg`)
- List of geometric nodes (corners)

Tips:
- In the header file, write the function `doc_structured_mesh` which creates these arrays.
- There are several ways (for instance squares can be split SW to NE or SE to NW) a structured mesh can be built, any will do.

Testing:
- In the unit test, call the function then feed the outputs to a [MetrisAPI](../../../src/API/MetrisAPI.hxx) object, and then that to a [MetrisRunner](../../../src/MetrisRunner/MetrisRunner.hxx) object, to get a [MeshBase](../../../src/Mesh/MeshBase.hxx) out of it:

```
MetrisAPI myAPI(...); // see MetrisAPI.hxx
MetrisParameters param; // can be left default
MetrisRunner(myAPI, param);
MeshBase &myMsh = *(run.msh_g);
writeMesh("name", myMsh); // visualize the mesh e.g. using Vizir
```
- Test the neighbor arrays by comparing with `myMsh.fac2fac`, entries should all be identical as the mesh elements are not reordered in this initialization.


Notes:
- The function is implemented in a header file so it can later be used in other exercises, using the the `MetrisAPI` into `MetrisRunner`.
- Doing this in 3D is more difficult: the main idea is to split a grid of cubes, in a way that their boundary splits correspond between neighbor cubes. For instance, always split `x` aligned cube faces from minimum `y` and `z` to maximum `y` and `z`. Then, the cubes still need to be covered in tetrahedra using those face splits.