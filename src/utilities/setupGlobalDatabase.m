%-------------------------------------------------------------------------%
% globDat = setupGlobalDatabase(K,Mesh,ProblemData,RVEdata) creates a
% global database 'globDat' that stores primary and derived variables
% useful for post-processing, and for reverting to earlier converged step
% if current step fails to converge. The enables one to resume the
% time-stepping with reduction in step-sizes.
%
% INPUT:  K              -> Macro-scale stiffness matrix
%         Mesh           -> Macro-scale mesh
%         ProblemData    -> Macro-scale problem data
%         RVEData        -> RVE data
%
% OUTPUT: globDat        -> MATLAB struct containing primary and derived
%                           variables.
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   06.01.2021
%
%-------------------------------------------------------------------------%

function globDat = setupGlobalDatabase(K,Mesh,ProblemData,RVEdata)

% Some useful variables
Nnodes = Mesh.noNodes;
Nelems = Mesh.noElements;
Ndofs  = ProblemData.problem.dof;
Nsteps = ProblemData.solver.nsteps;

%%%%%%%%%%%%%%%% Post-processing Variables %%%%%%%%%%%%%%%%

% Store macro-scale force-displacement
globDat.Force = zeros(Nnodes * Ndofs,Nsteps+1);
globDat.Disp  = zeros(Nnodes * Ndofs,Nsteps+1);

% Homogenised dual quantities for each integration point in a macro-scale
% element. Note that one integration point is used, so total number of
% integration points = number of macro-scale elements.
% (sig - stress in Voigt, dam - damage)
globDat.stress = cell(1,Nelems);
globDat.damage = cell(1,Nelems);


%%%%%%%%%%%%%%%% Previous Converged Step %%%%%%%%%%%%%%%%
% [ ** Serves as the restart point if current step do not converge]

% Old Step Stiffness Matrix
globDat.oldStepK         = K;

% Old Step Residuals
globDat.oldStepRes       = zeros(Nnodes * Ndofs,1);

% Old Step Solution
globDat.oldStepSol       = zeros(Nnodes * Ndofs,1);

% Old Step RVE Solution
globDat.RVEoldStepSol    = cell(1,Nelems);
for i = 1:Nelems
    globDat.RVEoldStepSol{i} = zeros(RVEdata.ndofs,1);
end

end

