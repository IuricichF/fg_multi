# FG_Multi

This library can be used for computing a discrete vector field preserving the multiparameter persistent homology of a multiparameter filtration.

The expected input is an .off file representing a simplicial complex (up to dimension three) with multiple scalar values defined at its vertices.

Example `mesh.off` (if you use this remove the commented lines form the input file)
```python
OFF             #keyword
8 12 0          #number of vertices - number of top simplices
0 1 0 1 2       #0 1 0 are the (x,y,z) coordinates, 1 2 are the vertex scalar values
1 1 -1 4 5
1 1 1 3 5
2 0 0 3 1
2 2 0 2 1
3 1 -1 6 0
3 1 1 9 9
4 1 0 1 2
3 0 1 3         #list of triangles 3 is the number of vertices 0 1 3 are the indexes of the correspinding vertices composing the triangle
3 0 2 3
3 0 1 4
3 0 2 4
3 1 3 5
3 2 3 6
3 1 4 5
3 2 4 6
3 3 5 7
3 3 6 7
3 4 5 7
3 4 6 7
```

If no scalar value is defined of the vertices the program will ask it should use the z coordinates or the x and y coordinates combined. Since the program can be used with an arbitrary number of scalar functions (even 1) it can be used for computing the Forman gradient of a filtration.

The output will be two txt files encoding the boundary matrices of the original simplicial complex and of the computed gradient.

Example `boundaryForman.txt`

```python
2                         #number of filtrations
1 2 3 6 9                 #scalar values (grades) for the first filtration
0 1 2 5 9                 #scalar values (grades) for the second filtration
0 : 0 2                   #0 dimension of the simplex, 0 2 indices of the filtration values in the list of grades
0 : 0 2
0 : 1 1
0 : 2 1
0 : 3 0
1 0 2 : 1 2               #1 dimension of the simplex, 0 2 indices of the boundary simplices, 1 2 indices of the filtration values in the list of grades
1 1 2 : 1 2
1 0 3 : 2 2
1 1 3 : 2 2
1 2 4 : 3 1
1 3 4 : 3 1
2 6 8 9 10 : 3 2
2 5 7 9 10 : 3 3
2 5 6 7 8 : 4 4
```

### Usage and Requirements

`Boost` is the only required library

From the main folder type the following command to compile the library

```bash
mkdir build
cd build
cmake ..
make
```

At this point, given an input file `mesh.off` you can run the program as follows

```bash
./NM_Forman mesh.off
```
