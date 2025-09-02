[Back to main.md](../main.md)

# Surface topology check

This function is implemented in [bunit/doc_surface_check.cxx](../../../bunit/doc_surface_check.cxx).

Inputs:
- Triangle and edge connectivity arrays `edg2poi`, `fac2poi`
- Neighbor arrays `edg2edg`, `fac2fac` (possibly `edg2fac`). 
- Element reference arrays `edg2ref`, `fac2ref`
- Boundary point information arrays `poi2bpo`, `bpo2ibi`
- A `MeshBase` object to call boundary edge hash table function `getedgglo()` (in `aux_topo.hxx`)

Outputs:
- Return whether the surface refs are correct:
  - Every boundary entity (triangles and edges in 3D) has a >= 0 reference
  - Two neighbor boundary entities have different references iff there lies a lower-dimensional entity inbetween (edge between triangles, node/corner between edges)
- Return whether point boundary information is correct:
  - Every point on the boundary has >= 0 `poi2bpo` entry
  - For `ipoin` on the surface, `ibpoi = poi2bpo[ipoin]` is the entry point to the point's boundary information in `bpo2ibi`. This is a linked list with entries `bpo2ibi(ibpoi,:) = ipoin, tdim, ientt, ibnxt` with `tdim` is the dimension (0, 1 or 2) of an entity `ientt`, and `ibnxt` is the next entry in `bpo2ibi` for ipoin. When there is no next entry, `ibnxt < 0`. 
  - Entries in `bpo2ibi` are sorted by topological dimension. 
  - If a point is dimension `pdim` (lowest dimensional entity attached), it lists all dimension `tdim > pdim` entities in its bpo2ibi entries. It lists only one entity of dimension `tdim = pdim`.


Notes:
- Edges can have >= 2 neighbors to a vertex. In that case, they store `-ienei - 1` in range `-infty, -2`. This forms a linked list (the neighbor then points to the next and so on, back to the original edge). This only happens if the vertex is a node/corner, as a curve can't meet itself in this way. 
- To know if a pair of vertices forms an edge in the mesh, call `getedgglo()`. This is constant time (if the hash table is well balanced) to query and insert into. 
- There is no global list of corners, edge neighbor information or `bpo2ibi` can be used to detect one. 

Testing:
- The structured mesh built in the [structured mesh exercise](structured_mesh.md) can be used to test the function. 
- A mesh can be read and some references made invalid, or entries in `bpo2ibi`

