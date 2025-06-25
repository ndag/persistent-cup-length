
## Example 1: Torus

We begin by showing how to use these functions on a simplicial model of a torus. 

First we specify the simplices needed for this model. These are written as objects of the class oin.Simplex which requires the specification of the vertices (listed in increasing order) and the time at which the simplex is added in the filtration. We put these simplices in a list and create a filtration from oin.Filtration(Simplices).

These simplex objects in the filtration also have the attribute id and the $0$-simplices should be listed so that their vertex matches their simplex id.

```python
#simplicial model of a torus.

#vertices
v0 = oin.Simplex([0],0.1)
v1 = oin.Simplex([1],0.2)
v2 = oin.Simplex([2],0.3)
v3 = oin.Simplex([3],0.4)
v4 = oin.Simplex([4],0.5)
v5 = oin.Simplex([5],0.6)
v6 = oin.Simplex([6],0.7)
v7 = oin.Simplex([7],0.8)
v8 = oin.Simplex([8],0.9)

#edges
e1 = oin.Simplex([0,1],1.1)
e2 = oin.Simplex([0,2],1.1)
e3 = oin.Simplex([0,3],1.1)
e4 = oin.Simplex([0,4],1.2)
e5 = oin.Simplex([0,6],1.2)
e6 = oin.Simplex([0,8],1.2)
e7 = oin.Simplex([1,2],1.3)
e8 = oin.Simplex([1,4],1.3)
e9 = oin.Simplex([1,5],1.3)
e10 = oin.Simplex([1,6],1.4)
e11 = oin.Simplex([1,7],1.4)
e12 = oin.Simplex([2,3],1.4)
e13 = oin.Simplex([2,5],1.5)
e14 = oin.Simplex([2,7],1.5)
e15 = oin.Simplex([2,8],1.5)
e16 = oin.Simplex([3,4],1.6)
e17 = oin.Simplex([3,5],1.6)
e18 = oin.Simplex([3,6],1.6)
e19 = oin.Simplex([3,7],1.7)
e20 = oin.Simplex([4,5],1.7)
e21 = oin.Simplex([4,7],1.7)
e22 = oin.Simplex([4,8],1.8)
e23 = oin.Simplex([5,6],1.8)
e24 = oin.Simplex([5,8],1.8)
e25 = oin.Simplex([6,7],1.9)
e26 = oin.Simplex([6,8],1.9)
e27 = oin.Simplex([7,8],1.9)

#triangles

t1 = oin.Simplex([0,1,4], 2.2)
t2 = oin.Simplex([0,1,6], 2.2)
t3 = oin.Simplex([0,2,3], 2.2)
t4 = oin.Simplex([0,2,8], 2.4)
t5 = oin.Simplex([0,3,4], 2.4)
t6 = oin.Simplex([0,6,8], 2.4)
t7 = oin.Simplex([1,2,5], 2.6)
t8 = oin.Simplex([1,2,7], 2.6)
t9 = oin.Simplex([1,4,5], 2.6)
t10 = oin.Simplex([1,6,7], 2.8)
t11 = oin.Simplex([2,3,5], 2.8)
t12 = oin.Simplex([2,7,8], 2.8)
t13 = oin.Simplex([3,4,7], 2.8)
t14 = oin.Simplex([3,5,6], 2.8)
t15 = oin.Simplex([3,6,7], 2.9)
t16 = oin.Simplex([4,5,8], 2.9)
t17 = oin.Simplex([4,7,8], 2.9)
t18 = oin.Simplex([5,6,8], 2.9)

Simplices = [v0, v1, v2, v3, v4, v5, v6, v7, v8, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19, e20, e21, e22, e23, e24, e25, e26, e27, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18]

fil = oin.Filtration(Simplices)
```

Now that we have a filtration we use the functions in Oineus to get the VRU Decomposition of this filtration. We will need the V and R matrices of this filtration which are written in a sparse matrix notation ([[0], [1], [0,2]] encodes the matrix with 1s in the [1,1], [2,2], [3,3] and [1,3] coordinates and zeros elsewhere). For our purposes, these matrices have coordinates in $\mathbb{Z}_2$

We note that we will want to use the dualized form of this to get Cocycle Representatives which reverses the order in which our cosimplices are listed i.e. if we have $n$ Simplices in Total, the dual simplex to simplex $i$, is indexed as $n-1-i$ in the matrices V and R. To refer to these cosimplices we create a list sorting our simplices in order of decreasing id.

```python
# cohomology
dualize = True
# create VRU decomposition object, does not perform reduction yet
dcmp = oin.Decomposition(fil, dualize)

rp = oin.ReductionParams()

rp.compute_u = rp.compute_v = True
rp.n_threads = 1
# perform reduction
dcmp.reduce(rp)

# now we can acess V, R and U
# indices are sorted_ids of simplices == indices in fil.cells()
V = dcmp.v_data
R = dcmp.r_data

simplices = fil.simplices()
cosimplices = sorted(simplices, key= lambda x: x.id , reverse = True)

```

### Persistent Cup-Length

We now go over how to use this information to calculate the Persistent Cup-Length Diagram associated to the cocycle representatives gained from this R and V Matrix.

Following Algorithm 2 in "Persistent Cup-Length" we need A dimension bound, which is implicit in explicit filtrations, the ordered list of cosimplices (above) the column reduced R Matrix (above) and the barcodes annotated by their representative cocycles listed first in increasing order of death times, then in the increasing order of birth times. 

We get these annotated barcodes by the get_Cocycle_Reps function which pulls out all of the barcodes from the oineus decomposition and pulls out the corresponding column of the V matrix which lists the cosimplices that compose the cocycle representing this barcode. To see how to get the persistent cohomology diagrams using Oineus go to https://github.com/anigmetov/oineus/blob/master/doc/tutorial.md

```python
RCocycles = get_Cocycle_Reps(dcmp,fil,V,2)
#We need to reorder these first by increasing death time then by increasing birth time
ReorderedCocycles = sorted(RCocycles, key= lambda x: (x.death, x.birth), reverse=False)
```

Finally we use these ordered cocycle representatives to get the persistent cup length diagram in the form of a Matrix whose rows indicate the birth times, the columns indicate the death times and the value in the i,j entry represents the value that should be assigned to that interval in the persistent cup length diagram

```python
Cup_diagram = Persistent_Cup_Diagram(2, cosimplices, R, ReorderedCocycles)

Cup_diagram_Matrix = Cup_diagram[0]
b_time = Cup_diagram[1]
d_time = Cup_diagram[2]
```
