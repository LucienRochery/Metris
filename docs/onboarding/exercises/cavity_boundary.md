[Back to main.md](../main.md)

# Cavity boundary function

This function is implemented in [bunit/doc_cavity_boundary.cxx](../../../bunit/doc_cavity_boundary.cxx). 

Inputs:
- A list of elements `lcavel` (or some other name) of topological dimension `tdim`.
- Element connectivity array: `fac2poi` in Metris
- Element neighbours array: `fac2fac` in Metris
- A tag `itag` and element array `ltag` such that, on function entry, `ltag(ielem) <= tag` for all `ielem`: `tag[ithread]` and `fac2tag(ithread,:)` in Metris

Outputs: 
- A list of unique facets on the boundary of the cavity: an `intAr2` object sized some `N x tdim`, with each entry being `tdim` vertices that comprise the facet of only one cavity element. 

Tips: 
- See [ball.md](ball.md) for tagging. 
- The main idea is to colour cavity elements first, then use that to check if a facet is shared between two cavity elements or belongs to only one. 


Notes:
- Meshes can be assumed manifold, i.e. a triangle only has one neighbour per edge. 
- A more sophisticated version can be implemented which takes a MshCavity object and builds a boundary list for each topological dimension (1, 2 and 3). 