[Back to main.md](../main.md)

# Cavity boundary function

A **cavity** is simply a list of elements, which can be of different topological dimension (although this exercise assumes single dimension). 
The **cavity operator** is the single operator for topological change in Metris, and works by reconnecting a cavity's boundary to a specified point. 
For this, it is necessary to be able to gather the cavity's boundary. 


This function is implemented in [bunit/doc_cavity_boundary.cxx](../../../bunit/doc_cavity_boundary.cxx). 

Inputs:
- A list of elements `lcavel` (or some other name) of topological dimension `tdim`.
- Element connectivity array: `fac2poi` in Metris
- Element neighbors array: `fac2fac` in Metris
- A tag `itag` and element array `ltag` such that, on function entry, `ltag(ielem) <= tag` for all `ielem`: `tag[ithread]` and `fac2tag(ithread,:)` in Metris

Outputs: 
- A list of unique facets on the boundary of the cavity: an `intAr2` object sized some `N x tdim`, with each entry being `tdim` vertices that comprise the facet of only one cavity element. 

Tips: 
- See [ball.md](ball.md) for tagging. 
- The main idea is to color cavity elements first, then use that to check if a facet is shared between two cavity elements or belongs to only one. 
- Use `lnoed2`, `lnofa3`  to get vertices of a triangle edge, resp tetrahedron face. See `metris_constants.hxx`.


Notes:
- Meshes can be assumed manifold, i.e. a triangle only has one neighbor per edge. 
- A more sophisticated version can be implemented which takes a MshCavity object which holds an edge, triangle and tetrahedron cavity, and builds a boundary list for each topological dimension. 
- There is currently no Metris function that implements this functionality, instead the boundary of a cavity is often looped over directly, in the same way as the solution to this exercise does to fill the array. 