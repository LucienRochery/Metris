[Back to main.md](../main.md)

# Bézier evaluation

This function is implemented in [bunit/doc_localization.cxx](../../../bunit/doc_localization.cxx).

Inputs:
- Control points of an element `ielem` of degree `ideg` and dimension `tdim` 
- Barycentric coordinates: `tdim + 1` doubles.

Outputs:
- Physical coordinates: `tdim` doubles.

We'll implement two different algorithms. 
Firstly, a generic evaluation algorithm that computes the Bernstein functions at the barycentric coordinate and multiplies against the degrees of freedom (control points). 
Secondly, the de Casteljau specialized recursive algorithm. 

Notes:
- High-order constants can be found in `ho_constants.hxx`: `getnnode(tdim, ideg)` returns number of nodes of element type, `ordfac.s[ideg][inode][:]` holds the multi-index of node `inode` (e.g. vertices are `ideg, 0, 0`, `0, ideg, 0`,  edge control points are `i1, 0, i2` with `i1+i2=ideg`, etc.), `mul2nod(i1,i2,...)` functions are the inverse of `ord(edg/fac/tet)` i.e. they return `inode` from a multi-index. 
- Other ordering arrays are found in `metris_constants.hxx` such as `lnoed2/3` and `lnofa3` which give the vertices associated to edges or faces. 
- These functions can be specialized to a given dimension (e.g. edges, triangles or tetrahedra). 

Testing:
- Both versions below can be implemented and the results compared. 
- They can also be benchmarked by using `get_wall_time()` 

## Generic evaluation

Given barycentrics `bary`, the physical coordinate is: 

`coop = sum_inode Bernstein(inode, bary) * control_point(inode)`

Hence we first implement functions `Bernstein(inode, bary)` (or better yet, `Bernstein(i1, i2, ..., bary)` taking multi-indices): 

`Bernstein(i1, i2, ..., bary) = (i1 + i2 ...)!  / (i1! i2! ... ) bary[0]^i1 * ...`

## de Casteljau algorithm

### Recursive relation (Bernsteins)

The de Casteljau algorithm relies on the following recurrence relationship, here for 1D Bernsteins:

$B_{ij}(\xi) = \xi_1 B_{i-1,j}(\xi) + \xi_2 B_{i,j-1}(\xi)$

with the convention that $B_{ij} = 0$  if $i < 0$ or $j < 0$. 
Indeed, in the general case,

$ B^d_{\alpha}(\xi) 
= \frac{d!}{\alpha!} \prod_{i=1}^{k+1} \xi_i ^{\alpha_i} 
= \frac{d}{\alpha_j} \xi_j \frac{(d-1)!}{(\alpha - e_j)!} \prod_{i=1}^{k+1} \xi_i^{\alpha_i - \delta_{ij}} \delta_{\alpha_i > 0}
= \frac{d}{\alpha_j} \xi_j B^{d-1}_{\alpha - e_j}(\xi)$

where:

- $j$ is any index between $0$ and $k+1$.
- $\alpha$ is a $(k+1)$-tuple that sums to $d$ (the superscript $d$ over $B$ is only for clarity as $\alpha$ specifies $d$ implicitely), 
- $\alpha! = \alpha_1!...\alpha_{k+1}!$,
- $e_j$ is the $(k+1)$-tuple of zeroes except for $1$ at the $j$-th position,
- and $\delta_{ij} = 0$ or $1$ iff $i=j$. 

We multiply this relation by $\alpha_j/d$ and sum over $j$, and obtain:

$  \sum_{j=1}^{k+1} \frac{\alpha_j}{d} B^d_{\alpha}(\xi) = \sum_{j=1}^{k+1}  \xi_j B^{d-1}_{\alpha - e_j}(\xi)$

But since $\alpha$ sums to $d$, we have $\sum_{j=1}^{k+1} \frac{\alpha_j}{d} = 1$, yielding the announced relation. 


### Recursive relation (element mappings)

This means the element mapping $F_K$ can be written as, e.g. here for a degree $d$ edge,

$$F_K(\xi) = \sum_{i + j = d} B_{ij}(\xi) P_{i_1j}
= \sum_{i + j = d} (\xi_1 B_{i-1,j}(\xi) P_{ij} + \xi_2 B_{i,j-1}(\xi) P_{ij})$$

We don't want to evaluate Bernstein polynomials, but rather to arrive at degree $d-1$ element evaluations. 
The trick is now to split the sum and do a change of variables, for instance for the first term:

$$
\sum_{0 \leq i,j \leq d, i + j = d} \xi_1 B_{i-1,j}(\xi) P_{ij} = 
\sum_{0 \leq j \leq d, 0 \leq i' + 1\leq d,  i' + 1 + j = d} \xi_1  B_{i',j}(\xi) P_{i'+1j} \\
$$

Now, notice that if $i' = -1$, then $B_{i'j} = 0$ hence $0 \leq i' + 1\leq d$ becomes $0 \leq i' \leq d-1$. 
Similarly, the constraint $i' + 1 + j = d$ forbids $j = d$ thus $0 \leq j \leq d$ becomes $0 \leq j \leq d-1$. 
Thus:
$$
\sum_{0 \leq i,j \leq d, i + j = d} \xi_1 B_{i-1,j}(\xi) P_{ij} 
= \sum_{0 \leq i,j \leq d-1, i + j = d - 1} \xi_1 B_{i,j}(\xi) P_{i+1j} 
$$


We notice the right hand side is simply the mapping of a degree $d-1$ edge, which has control points taken from a subset of those of the original element. 

Algorithmically, the easiest way to carry this out is to  pass in an index offset (representing the +1 in i+1 above) to the eval function, which then calls itself for a lower degree in $k+1$ ways, each with +1 offset at a different index.  

For instance, for edges:

```
eval1(ideg, nodes, bary, di, dj){
  if(ideg > 1){
    return bary[0] * eval1(ideg - 1, nodes, bary, di + 1, dj    )
         + bary[1] * eval1(ideg - 1, nodes, bary, di    , dj + 1);
  }else{
    int i10 = mul2nod(1+di, 0+dj);
    int i01 = mul2nod(0+di, 1+dj);
    return bary[0] * nodes[i10] + bary[1] * nodes[i01];
  }
}
```
