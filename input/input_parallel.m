% input_parallel.m contains the necessary input for parallel execution of
% the RVE problems in COMSOL Multiphysics.
%
% EXAMPLE: If you would like to use 4 physical cores, you can run 2 RVEs in
%          parallel, alloting 2 cores for each RVE problem. This is
%          achieved using:
%          inp_parallel.nparRVE  = 2;
%          inp_parallel.ncoreRVE = 2;


% Number of RVEs to solve in parallel
inp_parallel.nparRVE = 4;

% Number of cores per RVE
inp_parallel.ncoreRVE = 1;