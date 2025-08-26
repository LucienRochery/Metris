[Back to main.md](../main.md)

# Neighbors construction

Neighbor (or adjacency) arrays link the elements in a mesh with each other, and are the basis for navigating a mesh. 
A naive implementation involves a nested loop over elements, which is quadratic in time: much too expensive. 
Using hash tables, this can be brought down to linear time and space complexity. 

This function is implemented in [bunit/doc_neighbors.cxx](../../../bunit/doc_neighbors.cxx).

Inputs:
- Element connectivity array: e.g. `fac2poi` in Metris

Outputs:
- Same size neighbors array e.g. `fac2fac` which, at `(iface,i)` stores the neighbor to `iface` opposite its `i`-th vertex, or `-1` if none. 

Tips:
- Two elements are neighbors if the facet they share is seen twice in the mesh. Looping over facets in the mesh, seek them in a facet hash table: if not in, then add with value the host element, otherwise take the previous host element and make it the current's neighbor. 
- `types.hxx` has some convenience typedefs for hash table types, `HshTab_I2I` and `HshTab_I3I` are the most appropriate: they have key a pair or triplet of integers and value a single integer. 
- Use `stup2()` and `stup3()` (`utils/aux_misc.hxx`) to get keys of the right C++ type from loose integers. These functions also order the values, which avoids making distinct keys out of e.g. `(2, 3)` and `(3,2)`. 
- Use `lnoed2`, `lnofa3`  to get vertices of a triangle edge, resp tetrahedron face. See `metris_constants.hxx`.

Testing:
- The structured mesh built in `structured_mesh.md` can be used to validate by comparing neighbors arrays. Since the neighbors are ordered per element, as long as the elements are ordered the same, a simple equality check of all entries will do. 
- The neighbors can also be compared to the `MeshBase`'s `edg2edg`, `fac2fac`, `tet2tet`. 