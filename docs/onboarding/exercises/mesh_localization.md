[Back to main.md](../main.md)

# Mesh localization

Localizing a point in a mesh given its physical coordinates consists in finding the mesh element in which it lands (if any), as well as its barycentric coordinates in that element. 
In the context of Metris, this is most often used for metric interpolation: a new vertex is localized in the back mesh, and the barycentric coordinates are used to compute the metric. 

This requires completing [P1 element localization](P1_element_localization.md) first. 

This function is implemented in [bunit/doc_localization.cxx](../../../bunit/doc_localization.cxx).

Inputs:
- Coordinates `coop` of a point 
- Element connectivity array: `fac2poi` or `tet2poi` in Metris
- Element neighbors array: `fac2fac` etc. 
- Any seed element `iele0` to start the localization. 

Outputs:
- Index of mesh element `ielem` such that `coop` lies inside it, and its barycentric coordinates. 


Notes:
- Geometric and topological dimension are the same. That is, if `coop` has `n` entries, then we're considering topological dimension `n` elements. Otherwise this becomes a boundary projection function. 
- In practice, localization is not seeded from random elements, Metris keeps track of a `poi2bak` array to seed localizations from as close as possible to the final element. However if only convex domains are considered, this can be avoided (at the cost of CPU time). 

Tips:
- We know how to get the barycentric coordinates of `coop` in a single element from the [P1 element localization](P1_element_localization.md) exercise. Otherwise, function `inventP1()` in `low_localization.hxx` can be used. 
- If `coop` has negative `i`-th barycentric coordinate in an element `ielem`, then it lies in the half-space aligned with `ielem`'s `i`-th facet and opposite its `i`-th vertex. Thus, one moves closer to `coop` by going to the `i`-th neighbor of `ielem`. 
- In some cases (depending on mesh configuration), moving to neighbors in this manner can lead to infinite loops if elements are not tagged once visited. 
- There can be up to `tdim` negative barycentric coordinates, hence several possible neighbors to visit. 

Testing:
- The structured mesh built in the [structured mesh exercise](structured_mesh.md) can be used to test the function. 
- Non-convex domains shouldn't be attempted as the choice of the seed `iele0` becomes crucial. 

