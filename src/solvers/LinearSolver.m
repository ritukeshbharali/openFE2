%-------------------------------------------------------------------------%
% [x,data] = LinearSolver(A,b,solverinfo) solves a linear system of
% equations Ax = b, x being the unknown. The MATLAB struct 'solverinfo'
% must contain the type of solver 'linsolver' and additional information
% pertaining to preconditioners (for iterative solvers). 
%
% INPUT:  A              -> Matrix
%         b              -> Vector
%         solverinfo     -> MATLAB struct containing solver choices
%
% OUTPUT: x              -> Solution vector
%         data           -> Solver related data (useful for monitoring the
%                          performance of iterative solvers/preconditioners
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   08.01.2021
%-------------------------------------------------------------------------%

function [x,data] = LinearSolver(A,b,solverinfo)

% Solve the linear problem
switch solverinfo.linsolver
    
    % Direct solver (suitable for small problems)
    case 'direct'
        x    = A\b;
        data = [];
        
    case 'gmres-amg'
        % HSL-AMG preconditioning
        
        % GMRES solver
        
    otherwise
        error('Wrong Input!')
end

end

