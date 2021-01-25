%-------------------------------------------------------------------------%
% [K,D,RVEdata] = AssembleInitMacroSystem(Mesh,ProblemData,RVEdata,
% comsolPort) initialises an RVE for each macro-scale element (one-point
% Gauss integration scheme) and computes the initial macro-scale stiffness
% matrix 'K' using the homogenised tangent moduli 'D' obtain using
% sensitivity analysis.
%
% INPUT:  Mesh           -> Macro-scale mesh
%         ProblemData    -> MATLAB struct containing macro-scale problem
%                           information
%         RVEdata        -> MATLAB struct containing RVE problem
%                           information
%         comsolPort     -> Vector containing Port numbers associated with 
%                           the COMSOL servers
%
% OUTPUT: K              -> Macro-scale stiffness matrix
%         D              -> Macro-scale homogenised tangent moduli
%         RVEdata        -> The input 'RVEdata' is returned with
%                           the RVE number of degrees of freedom
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   06.01.2021
%
%-------------------------------------------------------------------------%

function [K,D,RVEdata] = AssembleInitMacroSystem(Mesh,...
                                 ProblemData,RVEdata,comsolPort)

% Define some useful variables                             
NN  = Mesh.noNodes * ProblemData.dof;
n   = ProblemData.nodes * ProblemData.dof;   
macroSpaceDim = ProblemData.SpaceDim;

% Compute Gauss point locations and weights
gp = IntegrationScheme(ProblemData.noInt, macroSpaceDim);

% Pre-defined variables to store the homogenised tangent moduli
D       = cell(1,Mesh.noElements);
RVEDofs = cell(1,Mesh.noElements);

parfor idx = 1:Mesh.noElements
    t = getCurrentTask();
    taskID=t.ID;
    try
        mphstart(comsolPort(taskID));
    catch
    end
    import('com.comsol.model.*');
    import('com.comsol.model.util.*');
    [D{1,idx},...
          RVEDofs{1,idx}] = RVEInit(idx, macroSpaceDim, RVEdata);
end

% Store RVE total dofs (works only if one RVE is allotted to all elements!)
RVEdata.ndofs = RVEDofs{1,1};

% Pre-allocation for global stiffness matrix
K = sparse([], [], [], NN, NN, 20*NN);

% Loop over elements
for idx = 1: Mesh.noElements
    
    % Set-up local node coordinates and connectivity
    connect_element = Mesh.connect(:,idx);
    x_element       = Mesh.x(:,connect_element);
    
    % Assemble matrix based on 1D,2D or 3D macro-scale model
    switch macroSpaceDim
        
        % 1D macro-scale model
        case 1
        
        uDofMap         = connect_element;
        
        % Loop over integration points
        for gauss_point = gp
            
            % Compute generic shape functions
            [N, dN, j] = ShapeFunction(x_element, gauss_point, ProblemData.type);

            % Compute shape functions for displacement and phase-field
            [~, ~, ~, Bs] = ShapeFunctionMultiDim(N, dN, ProblemData);

            % Reference to real area mapping
            dX = gauss_point(1) * j;
                        
            % Compute the jacobian (displacement and phase-field)
            KU_e = Bs' * D{1,idx} * Bs * dX * ProblemData.thickness{1,idx};
            
        end
        
        % 2D macro-scale model
        case 2
            
            uDofMap           = zeros(n,1);
            uDofMap(1:2:end)  = 2*connect_element-1;
            uDofMap(2:2:end)  = 2*connect_element;
            
            % Loop over integration points
            for gauss_point = gp
                
                % Compute generic shape functions
                [N, dN, j] = ShapeFunction(x_element, gauss_point, ProblemData.type);

                % Compute shape functions for displacement and phase-field
                [~, BU, ~, ~] = ShapeFunctionMultiDim(N, dN, ProblemData);

                % Reference to real area mapping
                dX = gauss_point(1) * j;

                % Compute the jacobian (displacement and phase-field)
                KU_e = BU' * D{1,idx} * BU * dX * ProblemData.thickness{1,idx};
                
            end
        
        % 3D macro-scale model    
        case 3
            error('3D not yet implemented. Choose 1D or 2D!')
        otherwise
                error('Wrong Input! Choose 1,2 or 3.')
    end
    
    % Assemble the global stiffness matrix
    K (uDofMap,uDofMap)   = K (uDofMap,uDofMap)   + KU_e;
    
end

% Check dimensions of the assembled macro-scale stiffness matrix
assert(length(K(:,1))  == Mesh.noNodes * ProblemData.dof, 'Incorrect K dimensions!');
assert(length(K(1,:))  == Mesh.noNodes * ProblemData.dof, 'Incorrect K dimensions!');
    
disp('     -> Computed the initial stiffness matrix')
disp(['        Size: ',num2str(length(K(:,1))),'x',num2str(length(K(1,:)))])
                             
end

