[Back to main](main.md)

# Table of Contents

- [Raw Containers](#raw-containers)
  - [MeshArray1D\<T,INT1\>](#mesharray1dtint1)
  - [MeshArray2D\<T,INT1,INT2\>](#mesharray2dtint1int2)
  - [Relevant specific MeshArrays](#relevant-specific-mesharrays)
    - [intAr1](#intar1--mesharray1dintmetris_int1)
    - [intAr2](#intar2--mesharray2dintmetris_int1metris_int2)
    - [dblAr1](#dblar1--mesharray1ddoublemetris_int1)
    - [dblAr2](#dblar2--mesharray2ddoublemetris_int1metris_int2)

- [Containers in Action: MeshBase Attributes](#containers-in-action-meshbase-attributes)
  - [Connectivity arrays](#connectivity-arrays)
  - [Vertices](#vertices)
  - [Neighbour arrays](#neighbour-arrays)
  - [Boundary information](#boundary-information)
  - [Vertex boundary links](#vertex-boundary-links)
  - [Utilities](#utilities)
  - [Hash tables](#hash-tables)
  - [Initialization](#initialization)
  - [Entity creation / deletion](#entity-creation--deletion)
  - [Misc geometry / mesh info](#misc-geometry--mesh-info)
  - [Work/Tag tracking (internal)](#worktag-tracking-internal-use)
  - [Debug / Testing](#debug--testing)

---

# Raw Containers

Metris implements various containers, essentially n-D arrays, customized to be better suited or more efficient for the intendend applications than their counterparts in the standard library. **Implementation:** [aux_msharrays.hxx](../../src/Arrays/aux_msharrays.hxx) and related files.

There exists an important aspect regarding the size of these containers. On the one hand, we have the size of the allocated memory for the array; we refer to this as its **buffer**. On the other hand, we have the **logical size** which refers to the subset of the buffer that is considered used. Of course, the logical size is always less or equal than the size of the buffer.

**Note:** In what follows, method signatures are in *pseudo-code*. Arguments' types hints are only included when particularly useful.

## `MeshArray1D<T,INT1>`

[Back to Index](#table-of-contents)

1D container for elements of type `T`, indexed with integer of type `INT1` (e.g. `int32_t`, `int64_t`, etc).
Internally backed by a `std::shared_ptr<T[]>` with lightweight semantics (soft copies share memory).

### Construction
- `MeshArray1D()` -> empty array, no memory allocation.
- `MeshArray1D(n)` -> allocate `n` elements (i.e. buffer of size `n`), but also sets logical size to `n`.
- `MeshArray1D(n, T* a)` -> wraps existing user-managed buffer. No `std::shared_ptr<T[]>` backup.
- `MeshArray1D({list})` -> from initializer list. Both buffer and logical size of the same size as the list.
- Copy / move constructors supported.
- Can also be built from a `MeshArray2D` (flatten).

After construction point, let now `arr` be of type `MeshArray1D<T,INT1>`.

### Access
- `arr[i]` -> element access (with bounds checks in debug).
- Iterators: `for (T& x : arr) ...` supported.
- `arr.size()` -> size of buffer.
- `arr.get_n()` -> logical size.
- `arr.n1_` -> reference alias to current logical size.

### Capacity and resizing
- `arr.allocate(m)` -> grow underlying buffer to at least size `m`.
- `arr.set_n(n)` -> set current logical size. Grow buffer if requested size bigger than current buffer.
- `arr.inc_n()` -> increase logical size by 1, growing buffer if needed.
- `arr.free()` -> release memory and reset size to 0 (buffer and logical size).

### Modification
- `arr.fill(T x)` -> fill logical size with `x`.
- `arr.copyTo(MeshArray1D<T,INT1>&out, ncopy)` -> deep copy to `out` up to `ncopy` elements.
- `arr.stack(T val)` -> append `val` **at the end of logical size** and increment logical size by one, growing buffer if needed.
- `arr.pop()` -> decrement logical size by one and return what was the last element of logical size. No memory operations are performed; technically the element is still in the buffer but is not considered as part of the logical size.

### Ownership
- `arr.set_sp(n, shared_ptr<T[]> a)` -> reset buffer to external shared storage `a`.
- `arr.get_sp()` -> access internal shared pointer.
- Copy assignment is **soft copy** (shares storage).
- Move assignment transfers ownership.

### I/O
- `print(FILE*)` or `print(std::ostream&)` -> pretty-print contents.

[Back to Index](#table-of-contents)

## `MeshArray2D<T,INT1,INT2>`

[Back to Index](#table-of-contents)

2D container for elements of type `T`, indexed by integers of types `INT1` (rows) and `INT2` (columns/stride).
Internally backed by a `std::shared_ptr<T[]>` with soft-copy semantics (copies share memory).
Implements row-major layout with a configurable stride (number of columns in memory).

### Construction
- `MeshArray2D()` -> empty array, no memory allocation.
- `MeshArray2D(m, s)` -> allocate buffer of `m × s` elements. Logical size set to `m` rows.
- `MeshArray2D(n, s, T* a)` -> wraps existing user-managed buffer (unmanaged).
- Copy / move constructors supported.

### Access
- `arr[i]` -> returns pointer to row `i`. Thus, `arr[i][j]` gives element `(i,j)`.
- `arr(i,j)` -> access element `(i,j)` with bounds checks.
- `arr.size()` -> total buffer size = `m1 × stride`.
- `arr.size1()` -> number of rows in buffer.
- `arr.size2()` -> stride (columns in memory).
- `arr.get_n()` -> logical number of rows.

### Capacity and resizing
- `arr.allocate(m, s)` -> grow/reallocate buffer to `m × s`.
  If `s` changes, data is copied/truncated as needed.
- `arr.set_n(n)` -> set logical number of rows. Will reallocate if needed.
- `arr.inc_n()` -> increment logical row count by 1, growing buffer if needed.
- `arr.free()` -> release memory and reset sizes.

### Modification
- `arr.fill(T x)` -> fill entire buffer with `x`.
- `arr.fill(n, s, T x)` -> fill only first `n × s` entries with `x`.
- `arr.copyTo(MeshArray2D<T,INT1,INT2>&out, ncopy)` -> deep copy to `out` (up to `ncopy` rows).

### Ownership
- `arr.get_sp()` -> access internal shared pointer.
- Copy assignment is **soft copy** (shares memory).
- Move assignment transfers ownership.

### I/O
- `print(FILE*)` or `print(std::ostream&)` -> pretty-print contents in nested `[row: [...]]` format.

[Back to Index](#table-of-contents)

## Relevant specific MeshArrays

[Back to Index](#table-of-contents)

There are a handful important types arising from `MeshArray` and specific types for their template arguments. These are arrays holding integers; extensively used for connectivity data structures, and arrays holding real numbers; useful for coordinates.

#### `intAr1 = MeshArray1D<int,METRIS_INT1>`
#### `intAr2 = MeshArray2D<int,METRIS_INT1,METRIS_INT2>`
#### `dblAr1 = MeshArray1D<double,METRIS_INT1>`
#### `dblAr2 = MeshArray2D<double,METRIS_INT1,METRIS_INT2>`

For more, see [types_arrays.hxx](../../src/types_arrays.hxx).

> **Note:** `METRIS_INT1-2` are internally defined (typically `int32_t` or `int64_t`) and subject to change. See [types_arrays.hxx](../../src/types_arrays.hxx).

[Back to Index](#table-of-contents)

# Containers in Action: `MeshBase` Attributes
[Back to Index](#table-of-contents)

The most important type in Metris is `MeshBase` (found in [MeshBase.hxx](../../src/Mesh/MeshBase.hxx)), which is the fundamental representation of a mesh. Naturally, the most relevant data structures are members of this type.

## `MeshBase`

### Connectivity arrays
- `edg2poi, fac2poi, tet2poi`
- `curdeg, strdeg` (mesh degree)
- Entity counts: `nedge, nface, nelem`
- Generic helpers: `nentt`, `ent2poi`

[Back to Index](#table-of-contents)

### Vertices
- `npoin`, `coord`
- `idim` (geometric dimension)
- `poi2ent`
- `poicstr`, `bb` (bounding box)
- Helpers: `getveredg`, `getverfac`, `getvertet`, `getverent`

[Back to Index](#table-of-contents)

### Neighbour arrays
- `edg2edg, fac2fac, fac2tet, tet2tet`
- `edg2fac`
- Generic helpers: `ent2ent`
- Lower-dim helpers: `fac2edg`, `tet2fac`

[Back to Index](#table-of-contents)

### Boundary information
- `CAD` (CADInfo)
- `edg2ref, fac2ref, tet2ref`
- `ndomn`, `is_nonmanifold`
- Helpers: `ent2ref`, `get_algnd`

[Back to Index](#table-of-contents)

### Vertex boundary links
- `poi2bpo`
- `bpo2ibi`
- `bpo2rbi`
- `nbpoi`
- Helper: `poi2ebp`

[Back to Index](#table-of-contents)

### Utilities
- Work arrays: `get_work`, `get_iwork`, `get_rwork`
- Tag arrays: `get_tagarray`, `ent2tag`, `ref2tag`
- Internal tags: `tag, tet2tag, fac2tag, ...`

[Back to Index](#table-of-contents)

### Hash tables
- `edgHshTab`, `facHshTab`

[Back to Index](#table-of-contents)

### Initialization
- `readConstants`, `copyConstants`
- `readMeshFile`, `readMeshData`
- `iniNeighbours`, `iniBdryPoints`, `iniCADLink`

[Back to Index](#table-of-contents)

### Entity creation / deletion
- `newpoitopo`, `newfactopo`, `newfacvirtual`, `newedgtopo`, `newedgvirtual`, `newbpotopo`
- `killpoint`, `rembpotag`

[Back to Index](#table-of-contents)

### Misc geometry / mesh info
- `facedg2glo`, `tetfac2glo`, `tetedg2glo`
- `get_geodev`
- `getBasis`, `setBasis`, `forceBasisFlag`
- `get_tdim`
- Periodicity: `isperiodic_face`, `nperiodic_face`

[Back to Index](#table-of-contents)

### Work/Tag tracking (internal use)
- iwork/rwork tracked arrays
- tagarrs, tagarr_locks
- debug helpers (`debug_get_*`)

[Back to Index](#table-of-contents)

### Debug / Testing
- `debug_get_iwork`, `debug_get_rwork`, `debug_get_work_lock`

[Back to Index](#table-of-contents)