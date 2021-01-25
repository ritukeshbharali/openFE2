# openFE2

openFE2 is a software framework capable of carrying out (linear/nonlinear) multi-scale FE2 simulations. The macro-scale problem is solved in MATLAB while the micro-scale (RVE) problems are solved in parallel in COMSOL Multiphysics. The parallel setup is achieved using the MATLAB Parallel Computing Toolbox, which creates multiple instances of the COMSOL Multiphysics server.

### Features
 - Beginner friendly functional programming architecture
 - Ability to restart from any saved state
 - Flexibility in adding new RVEs 
 - Flexibility in alloting physical cores for the RVE problems
 - GMSH v2, COMSOL or user-defined MATLAB macro-scale meshes
 
 ### Limitations
 - Tested only on Windows 10 operating system
 - Bar, triangle and tetrahedon elements for 1D/2D/3D macro-scale problems
 - Single Gauss point integration scheme on the macro-scale 

### Requirements
 - COMSOL Multiphysics 5.5
 - MATLAB R2019a and above, with Parallel Computing Toolbox
 
 ### Getting started
 - Install COMSOL Multiphysics 5.5 and MATLAB with the Parallel Computing Toolbox
 - Add path to the 'mli' library directory
 - Add path to the 'comsolmphserver.exe'
 - Provide suitable inputs in the input folder and run script 'openFE2.m'
 
 ### License
 [MIT](https://choosealicense.com/licenses/mit/)
