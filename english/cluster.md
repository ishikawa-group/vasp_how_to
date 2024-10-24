# Making nanocluster with ASE

## Making cluster
```python
from ase.cluster import Icosahedron
from ase.visualize import view
from ase.io import write

cluster = Icosahedron("Pt", noshells=5)
cluster.cell = [30, 30, 30]
cluster.center()

write("cluster.pdb", cluster)
```

### Loading cluster
```python
from ase.io import write, read

cluster = read("cluster.xyz")
cluster.cell = [30, 30, 30]
write("cluster.pdb", cluster)
```

## Distributing water around cluster (using packmol)
* Install packmol following the official website.

### Input files
* water pdb file: `water.pdb`
```
HEADER    water
COMPND
SOURCE
HETATM    1  H   HOH     1       9.626   6.787  12.673
HETATM    2  H   HOH     1       9.626   8.420  12.673
HETATM    3  O   HOH     1      10.203   7.604  12.673
CONECT    1    3
CONECT    2    3
CONECT    3    1    2
END
```

* packmol input file: `input.inp`
```
tolerance 2.3
filetype pdb
output output.pdb

structure cluster.pdb
  number 1
  fixed 0.0 0.0 0.0  0.0 0.0 0.0
end structure

structure water.pdb
  number 200
  inside box  2.0 2.0 2.0 28.0 28.0 28.0
end structure
```
* `inside box`: change according to your unit cell size of cluster

### Execute packmol
* `packmol < input.inp`
  + `output.pdb` is generated
