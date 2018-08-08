# meshgen

Generate mesh for 2D composite material using MSHR https://bitbucket.org/fenics-project/mshr

Meshes generations in composite.py [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/petch/meshgen/master?filepath=composite.ipynb)

Solve linear elasticity problem using FEniCS https://fenicsproject.org/

Three variants of solutions: 
* fine - simple FEM solution on fine mesh
* coarse - solution using numerical homogenization
* sparse - partial fine and coarse meshes

Elasticity solution in elasticity.py [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/petch/meshgen/master?filepath=elasticity.ipynb)

Tested on MSHR and FEniCS versions: 2017.2, 2018.1
