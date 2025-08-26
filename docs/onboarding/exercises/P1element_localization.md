[Back to main.md](../main.md)

# P1 element localization

Computing barycentric coordinates from physical coordinates (localization or inverse evaluation) is a basic geometric operation that is mainly used for localizing points in meshes (where the final element is initially unknown). 
In the case of linear elements, barycentric coordinates have a closed-form expression as a function of physical coordinates, so this is simple and very cheap to compute. 
This can be generalized to projection (e.g. where the point is 3D but the element is a triangle or edge) or to inverse evaluation in high-order elements (iterative process). 

This function is implemented in [bunit/doc_localization.hxx](../../../bunit/doc_localization.hxx) and tested in [bunit/doc_localization.cxx](../../../bunit/doc_localization.cxx).

Inputs:
- Coordinates `coop` of a point 
- An element `ielem`

Outputs:
- Barycentric coordinates `bary`, that is `tdim + 1` doubles (possibly outside the range [0,1]) of `coop` in `ielem`


Notes:
- Geometric and topological dimension are the same. For instance, if `coop` has 3 coordinates, then `ielem` is a tetrahedron (`tdim = 3`).

Tips:
- The `i`-th barycentric coordinate is the ratio of the measure (volume/area) of the element made by replacing the `i`-th vertex of `ielem` with `coop`, over the measure of `ielem`. 
- `linalg/det.hxx` has determinant functions (`detmat`) to help compute volumes

Testing:
- Functions `eval1`, `eval2` and `eval3` (`low_eval.hxx`) compute coordinates from barycentric coordinates. Use these to compare the result of `eval(bary)` with `coop`.
