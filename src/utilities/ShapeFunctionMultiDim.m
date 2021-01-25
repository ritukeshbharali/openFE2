%-------------------------------------------------------------------------%
% [Nv, Bv, Ns, Bs] = ShapeFunctionMultiDim(N, dN, ProblemData) returns the
% shape function and derivative for multi-dimensional system.
%
% INPUT:  N            - Shape Function
%         dN           - Derivative of Shape Function
%         ProblemData  - Problem Data struct, must contain 'SpaceDim',
%                       'nodes' and 'dof'
%
%
% OUTPUT: Nv          - Shape Function of a vector field
%         Bv          - Derivative of Shape Function of a vector field
%         Ns          - Shape Function of a scalar field
%         Bs          - Derivative of Shape Function of a scalar field
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   11.01.2021
%-------------------------------------------------------------------------%

function [Nv, Bv, Ns, Bs] = ShapeFunctionMultiDim(N, dN, ProblemData)

n = ProblemData.nodes;
m = ProblemData.dof;

if ProblemData.SpaceDim == 1 % 1D bar
    
    % Compute N-matrix for a scalar field
    Ns = N';
    
    % Compute the B-maxtrix for a scalar field
    Bs = dN';
    
    Nv = [];
    Bv = [];

elseif ProblemData.SpaceDim == 2
    
    % Compute the N-maxtrix for a vector field
    Nv = zeros(m,m*n);
    for i = 1:m
        Nv(i,i:m:end) = N(:)';
    end
    
    % Compute the B-maxtrix for a vector field
    Bv = zeros(3,m*n);
    for i = 1:m
        Bv(i,i:m:end) = dN(:,i);
    end
    Bv(3, 1:m:end) = dN(:,2)'; 
    Bv(3, 2:m:end) = dN(:,1)';
    
    % Compute N-matrix for a scalar field
    Ns = N';
    
    % Compute the B-maxtrix for a scalar field
    Bs = zeros(2,n);
    Bs(1,:) = dN(:,1);
    Bs(2,:) = dN(:,2);
    
elseif ProblemData.SpaceDim == 3
    error('3D not yet programmed!')
else
    error('Wrong input!')
end



end

