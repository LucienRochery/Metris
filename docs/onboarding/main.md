This file is a work in progress to help getting familiar with meshing and Metris. 

Directory exercises/ contains exercices to work through and solutions/ the corresponding descriptions of the algorithms/discussion. 

Related scripts are found in bunit/doc_*.cxx. They may eventually be moved to a subdirectory of docs/onboarding but that'll require a little bit of CMake to make convenient. 

Note: in all exercises, a `MeshBase& msh` object can be taken in instead of only the necessary loose arrays. 

# Topological functions

- [ball](exercises/ball.md): gather set of elements surrounding a vertex. 
- [cavity boundary](exercises/cavity_boundary.md): gather facets of cavity (list of elements) boundary.
- [structured mesh](exercises/structured_mesh.md): generate a structured mesh 
- [neighbors construction](exercises/neighbors_construction.md): generate the neighbors arrays using hash tables

# Geometry functions

- [P1 element localization](exercises/P1_element_localization.md): get the barycentric coordinates of a point in an element

# Classic algorithms

- [mesh localization](exercises/mesh_localization.md): find which mesh element a point (not vertex) lies in. 

# High-order meshes

- [degree elevation](exercises/degree_elevation.md): add degrees of freedom to a degree 1 mesh to make it a higher degree. 