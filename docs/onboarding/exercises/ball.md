[Back to main.md](../main.md)

# Ball function

This function is implemented and tested against the functions in `src/low_topo.hxx` in [bunit/doc_ball.cxx](../../../bunit/doc_ball.cxx). 

Inputs:
- A vertex `ipoin`
- Element connectivity array: `fac2poi` in Metris
- Element neighbours array: `fac2fac` in Metris
- A tag `itag` and element array `ltag` such that, on function entry, `ltag(ielem) <= tag` for all `ielem`: `tag[ithread]` and `fac2tag(ithread,:)` in Metris
- A seed element `iele0` which contains vertex `ipoin`. 

Outputs: 
- A list of unique entries, all the elements that contain vertex ipoin: an `intAr1` object in Metris (or `std::vector`).

Tips: 
- On entry, we will have `fac2tag(ithread,iface) == tag[ithread]` for some `iface` but not others: we must first increment `tag[ithread]`, and we can then do `fac2tag(ithread,iface) = tag[ithread]` to signify that an `iface` has been marked. 
- In the 2D case, and for a vertex in the interior of the domain, we can get the elements in a certain order (looping around the vertex). Implement this routine for the most general case where this is not true, as is the case in 3D or in 2D if the vertex is on the boundary. 
- The main idea is to seed the ball with `iele0` and then have a main loop going through the current ball and adding any neighbours which 1) contain the point and 2) have not yet been added to the ball. 


Notes:
- Dimension specific versions `doc_ball2()` and `doc_ball3()` can be implemented, and then a dimension generic one `doc_ball()` which gets all the elements of a provided topological dimension (use `msh.ent2poi`, `ent2tag` etc. helpers). 
- A function `doc_ball2_interior()` can be written which uses the fact triangles are ordered around a vertex. 
- Meshes can be assumed manifold, i.e. a triangle only has one neighbour per edge. 