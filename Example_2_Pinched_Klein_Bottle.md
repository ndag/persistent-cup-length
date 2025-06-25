This is the filtration given as Example 17 in "Persistent Cup-Length"

Start with the simplices that build the filtration

```python
#simplicial model of a Klein Bottle.

#vertices
w0 = oin.Simplex([0],0.1)
w1 = oin.Simplex([1],1.0)
w2 = oin.Simplex([2],1.0)
w3 = oin.Simplex([3],2.0)
w4 = oin.Simplex([4],2.0)
w5 = oin.Simplex([5],2.0)
w6 = oin.Simplex([6],2.0)
w7 = oin.Simplex([7],2.0)
w8 = oin.Simplex([8],2.0)

#edges
f1 = oin.Simplex([0,1],1.0)
f2 = oin.Simplex([0,2],1.0)
f3 = oin.Simplex([1,2],1.0)
f4 = oin.Simplex([0,3],2.0)
f5 = oin.Simplex([0,5],2.0)
f6 = oin.Simplex([0,6],2.0)
f7 = oin.Simplex([0,7],2.0)
f8 = oin.Simplex([1,3],2.0)
f9 = oin.Simplex([1,4],2.0)
f10 = oin.Simplex([1,5],2.0)
f11 = oin.Simplex([1,8],2.0)
f12 = oin.Simplex([2,3],2.0)
f13 = oin.Simplex([2,4],2.0)
f14 = oin.Simplex([2,7],2.0)
f15 = oin.Simplex([2,8],2.0)
f16 = oin.Simplex([3,4],2.0)
f17 = oin.Simplex([3,6],2.0)
f18 = oin.Simplex([3,8],2.0)
f19 = oin.Simplex([4,5],2.0)
f20 = oin.Simplex([4,6],2.0)
f21 = oin.Simplex([4,7],2.0)
f22 = oin.Simplex([5,6],2.0)
f23 = oin.Simplex([5,7],2.0)
f24 = oin.Simplex([5,8],2.0)
f25 = oin.Simplex([6,7],2.0)
f26 = oin.Simplex([6,8],2.0)
f27 = oin.Simplex([7,8],2.0)

#triangles

tr1 = oin.Simplex([0,1,2], 3.0)
tr2 = oin.Simplex([0,1,3], 2.0)
tr3 = oin.Simplex([0,1,5], 2.0)
tr4 = oin.Simplex([0,2,3], 2.0)
tr5 = oin.Simplex([0,2,7], 2.0)
tr6 = oin.Simplex([0,5,6], 2.0)
tr7 = oin.Simplex([0,6,7], 2.0)
tr8 = oin.Simplex([1,2,4], 2.0)
tr9 = oin.Simplex([1,2,8], 2.0)
tr10 = oin.Simplex([1,3,4], 2.0)
tr11 = oin.Simplex([1,5,8], 2.0)
tr12 = oin.Simplex([2,3,8], 2.0)
tr13 = oin.Simplex([2,4,7], 2.0)
tr14 = oin.Simplex([3,4,6], 2.0)
tr15 = oin.Simplex([3,6,8], 2.0)
tr16 = oin.Simplex([4,5,6], 2.0)
tr17 = oin.Simplex([4,5,7], 2.0)
tr18 = oin.Simplex([5,7,8], 2.0)
tr19 = oin.Simplex([6,7,8], 2.0)

Klein_Simplices = [w0, w1, w2, w3, w4, w5, w6, w7, w8, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22, f23, f24, f25, f26, f27, tr1, tr2, tr3, tr4, tr5, tr6, tr7, tr8, tr9, tr10, tr11, tr12, tr13, tr14, tr15, tr16, tr17, tr18, tr19]

Klein_fil = oin.Filtration(Klein_Simplices)

```
Then get the VRU Decomposition

```python
# cohomology
dualize = True
# create VRU decomposition object, does not perform reduction yet
Klein_dcmp = oin.Decomposition(Klein_fil, dualize)

rp = oin.ReductionParams()

rp.compute_u = rp.compute_v = True
rp.n_threads = 1
# perform reduction
Klein_dcmp.reduce(rp)

# now we can acess V, R and U
# indices are sorted_ids of simplices == indices in fil.cells()
Klein_V = Klein_dcmp.v_data
Klein_R = Klein_dcmp.r_data

Klein_simplices = Klein_fil.simplices()
Klein_cosimplices = sorted(Klein_simplices, key= lambda x: x.id , reverse = True)

```
Use the decomposition to pull out barcodes with their cocycle representatives and order them appropriately

```python
Klein_Rep_Cocycles = get_Cocycle_Reps(Klein_dcmp, Klein_fil, Klein_V, 2)

Klein_ReorderedCocycles = sorted(Klein_Rep_Cocycles, key= lambda x: (x.death, x.birth), reverse=False)

```

Then run these cocycle representatives through the Persistent Cup Length function 

```python
Klein_Cup_diagram = Persistent_Cup_Diagram(4, Klein_cosimplices, Klein_R, Klein_ReorderedCocycles)

Klein_Cup_diagram_Matrix = Klein_Cup_diagram[0]
Klein_b_time = Klein_Cup_diagram[1]
Klein_d_time = Klein_Cup_diagram[2]
```
