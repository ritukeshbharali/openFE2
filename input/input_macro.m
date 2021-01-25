% input_macro.m contains the necessary input for defining the macro-scale
% problem.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPATIAL DIMENSION, MESH, ELEMENT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spatial dimension
inp_macro.problem.SpaceDim = 1;

% Mesh
inp_macro.problem.mesh_file = 'MATLAB_mesh.m';

% Nodes per elements
inp_macro.problem.nodes = 2;

% DOFs per node
inp_macro.problem.dof = 1;

% Number of integration points per element
inp_macro.problem.noInt =  1;

% Element type (Available types: Bar2, Tri3, Tet4)
inp_macro.problem.type = 'Bar2';

% Element thickness
inp_macro.problem.thickness{1,1} = 0.9;
inp_macro.problem.thickness{1,2} = 1;
inp_macro.problem.thickness{1,3} = 1;
inp_macro.problem.thickness{1,4} = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRAINT DATA
%
% NOTE: Macro-scale problem is displacement-controlled. Therefore,
% 'constraints' must be an array (dofs,multiplier).
% 'multiplier' = 0 [homogeneous Dirichlet conditions]
%              = 1 [ihhomogeneous Dirichlet conditions]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inp_macro.constraints       = zeros(6,2);
% inp_macro.constraints(1,:)  = [1, 0];
% inp_macro.constraints(2,:)  = [2, 0];
% inp_macro.constraints(3,:)  = [3, 1];
% inp_macro.constraints(4,:)  = [5, 1];
% inp_macro.constraints(5,:)  = [7, 0];
% inp_macro.constraints(6,:)  = [8, 0];

inp_macro.constraints       = zeros(2,2);
inp_macro.constraints(1,:)  = [1, 0];
inp_macro.constraints(2,:)  = [5, 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVER DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of steps
inp_macro.solver.nsteps = 4;

% Initial step-size
inp_macro.solver.initstepsize = 5e-3;

% Terminating step-size
inp_macro.solver.finalstepsize = 1e-8;

% Cut step-size at these step number [must be an empty/filled array]
inp_macro.solver.cutsteps = [];

% Factor reduction in step-size [timestep = timestep/10 means factor 10]
inp_macro.solver.cutstepsize = 10;

% Fraction of peak load for termination
% Example: 0.25 means if current load is less than 25% of peak load,
% the simulation would terminate.
inp_macro.solver.loadfraction = 0.1;

% Non-linear Solver
inp_macro.solver.nonlinsolver = 'newton';

% Maximum number of iterations
inp_macro.solver.maxiter = 50;

% Stiffness Update
inp_macro.solver.stiffupdate = 'initial';

% Non-linear Solver Tolerance
inp_macro.solver.nonlinsolvertol = 1e-2;

% Linear Solver
inp_macro.solver.linsolver = 'direct';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESTART OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Steps at which restart files are stored [must be an empty/filled array]
inp_macro.restart.steps = [10 20 30 40];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Monitoring Dofs for load-displacement plot
inp_macro.postproc.loadDof = 1;
inp_macro.postproc.dispDof = 5;





