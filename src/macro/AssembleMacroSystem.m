%-------------------------------------------------------------------------%
% [K,Res,GPdata,Mflag] = AssembleMacroSystem(Mesh,Sol,MacroData,
%  RVEdata,comsolPort,step,iter) assembles the macro-scale residual and/or
%  the macro-scale stiffness matrix based user-defined choices pertaining
%  to the different variants of the Newton-Raphson method. This involves
%  solving Finite Element BVP at the RVE level for each macro-scale Gauss
%  point.
%
% INPUT:  Mesh           -> Macro-scale mesh
%         Sol            -> Macro-scale solution
%         MacroData      -> MATLAB struct containing macro-scale problem
%                           information
%         RVEdata        -> MATLAB struct containing RVE problem
%                           information
%         step           -> Current step number in the analysis
%         iter           -> Current iteration in the analysis
%         comsolPort     -> Vector containing Port numbers associated with 
%                           the COMSOL servers
%
% OUTPUT: K              -> Macro-scale stiffness matrix
%         Res            -> Macro-scale residual vector
%         GPdata         -> Macro-scale Gauss-point data (stresses)
%         Mflag          -> Flag indicating the success/failure (1/0) in
%                           computing the macro-scale stiffness matrix and
%                           residual. Structure [Stiffness Matrix, Residual]
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   08.01.2021
%
%-------------------------------------------------------------------------%

function [K,Res,GPdata,Mflag] = AssembleMacroSystem(...
                           Mesh,Sol,MacroData,RVEdata,comsolPort,...
                                               varargin)
                                           
% Define some useful variables
NN            = Mesh.noNodes * MacroData.problem.dof;
n             = MacroData.problem.nodes * MacroData.problem.dof; 
MacroSpaceDim = MacroData.problem.SpaceDim;

Mflag = [1,1];

% Compute Gauss point locations and weights
gp = IntegrationScheme(MacroData.problem.noInt, MacroSpaceDim);

% Check whether stiffness matrix assembly is required
switch nargin
    case 6
        action = varargin{1};
        switch action
            case 'Jac'
            case 'NoJac'
            otherwise
                error('Wrong input!')
        end        
    case 7
        step = varargin{1};        
        iter = varargin{2};  
        assert(isnumeric(step)==1,'Sixth input argument must be a number!');
        assert(isnumeric(iter)==1,'Seventh input argument must be a number!');
        
        switch MacroData.solver.stiffupdate
            
            % Initial stiffness method
            case 'initial'

                action = 'NoJac';

            % Modified Newton method    
            case 'step'

                if step == 2 && iter == 1
                    action = 'NoJac';
                elseif step > 2 && iter == 1
                    action = 'Jac';
                else
                    action = 'NoJac';
                end 

            % Full Newton method    
            case 'iter'

                if step == 2 && iter == 1
                    action = 'NoJac';
                else
                    action = 'Jac';
                end

            otherwise
                error('Wrong Input! Choose "initial", "step" or "iter"')
        end
            
        
    otherwise
        error('Unexpected number of input arguments!')
end

% Assemble (1D, 2D or 3D) macro-scale residual and/or stiffness matrix
switch action
    
    case 'NoJac'
        
        % COMPUTE MACRO-SCALE STRAINS       
        % Pre-allocation to store the macro-scale Gauss Point (GP) strains
        strain = cell(1,Mesh.noElements);
        
        % Compute the macro-scale GP strains
        for idx = 1:Mesh.noElements
            
           % Set-up local node coordinates and connectivity
            connect_element = Mesh.connect(:,idx);
            x_element       = Mesh.x(:,connect_element);

            % Get the local displacement and phase-field
            switch MacroSpaceDim
                
                case 1
                    uDofMap           = connect_element;                    
                    localU            = Sol(uDofMap);
                    
                    % Loop over integration points
                    for gauss_point = gp

                        % Compute generic shape functions
                        [N, dN, ~] = ShapeFunction(x_element, gauss_point, MacroData.problem.type);

                        % Compute shape functions for displacement and phase-field
                        [~, ~, ~, Bs] = ShapeFunctionMultiDim(N, dN, MacroData.problem);

                        % Get strain, pf and its gradient at integration point
                        % (gauss_point)
                        strain{1,idx} = Bs * localU;
                    end
                    
                case 2           
                    uDofMap           = zeros(n,1);
                    uDofMap(1:2:end)  = 2*connect_element-1;
                    uDofMap(2:2:end)  = 2*connect_element;
                    localU            = Sol(uDofMap);

                    % Loop over integration points
                    for gauss_point = gp

                        % Compute generic shape functions
                        [N, dN, ~] = ShapeFunction(x_element, gauss_point, MacroData.problem.type);

                        % Compute shape functions for displacement and phase-field
                        [~, BU, ~, ~] = ShapeFunctionMultiDim(N, dN, MacroData.problem);

                        % Get strain, pf and its gradient at integration point
                        % (gauss_point)
                        strain{1,idx} = BU * localU;

                    end
                case 3
                    error('3D not yet implemented!')
                otherwise
                    error('Wrong input!')
            end                   
        end
        
        % SOLVE THE RVE PROBLEMS
        % Pre-allocation to store the macro-scale Gauss Point (GP) stresses
        stress = cell(1,Mesh.noElements);
        
        parfor idx = 1:Mesh.noElements
            t = getCurrentTask();
            taskID=t.ID;
            try
                mphstart(comsolPort(taskID));
            catch
            end 
            import('com.comsol.model.*');
            import('com.comsol.model.util.*');
            
            [~,stress{1,idx},flag{1,idx}] = RVESolve(idx,strain{1,idx},...
                                         MacroSpaceDim,'NoJac',RVEdata);
        end
        
        % Store the homogenised stresses
        GPdata.stress = stress;
        
        
        % ASSEMBLE MACRO-SCALE RESIDUAL
        
        % Check if all RVE problems were successfully solved
        for idx = 1:Mesh.noElements
            el_flag = flag{1,idx};
            if el_flag(2) == 0
                res_flag = 0;
                break;
            else
                res_flag = 1;
            end
        end
        
        % Assemble only if all RVE problems were solved successfully
        if res_flag == 1
            
            % Pre-allocate the macro-scale residual vector
            Res = zeros(NN,1);
        
            % Loop over elements
            for idx = 1:Mesh.noElements

                % Set-up local node coordinates and connectivity
                connect_element = Mesh.connect(:,idx);
                x_element       = Mesh.x(:,connect_element);

                switch MacroSpaceDim

                    case 1 % 1D
                        uDofMap           = connect_element;

                        % Loop over integration points
                        for gauss_point = gp
                            % Compute generic shape functions
                            [N, dN, j] = ShapeFunction(x_element, gauss_point, MacroData.problem.type);

                            % Compute shape functions for displacement and phase-field
                            [~, ~, ~, Bs] = ShapeFunctionMultiDim(N, dN, MacroData.problem);

                            % Reference to real area mapping
                            dX = gauss_point(1) * j;

                            % Compute the residual (displacement)
                            Res_e = Bs' * stress{1,idx} * dX * MacroData.problem.thickness{1,idx};

                        end

                    case 2 % 2D               
                        % Get the local displacement and phase-field (2-D)
                        uDofMap           = zeros(n,1);
                        uDofMap(1:2:end)  = 2*connect_element-1;
                        uDofMap(2:2:end)  = 2*connect_element;

                        % Loop over integration points
                        for gauss_point = gp
                            % Compute generic shape functions
                            [N, dN, j] = ShapeFunction(x_element, gauss_point, MacroData.problem.type);

                            % Compute shape functions for displacement and phase-field
                            [~, BU, ~, ~] = ShapeFunctionMultiDim(N, dN, MacroData.problem);

                            % Reference to real area mapping
                            dX = gauss_point(1) * j;

                            % Compute the residual (displacement)
                            Res_e = BU' * stress{1,idx} * dX * MacroData.problem.thickness{1,idx};

                        end

                    case 3
                        error('3D not yet implemented!')
                    otherwise
                        error('Wrong input!')
                end                    

                % Assemble the global residuals
                Res (uDofMap)  = Res (uDofMap)  + Res_e;

            end

            % Check dimensions of the assembled macro-scale residual
            assert(length(Res)  == Mesh.noNodes * MacroData.problem.dof, 'Incorrect Res dimensions!');

            K       = [];
            Mflag(2) = 1;
            Mflag(1) = 0;
            
        elseif res_flag == 0
            
            Res = [];
            K   = [];
            Mflag(2) = 0;
            Mflag(1) = 0;
            
        else          
            error('Check .\src\macro\AssembleMacroSystem.m for errors!')
            
        end
            
        
    case 'Jac'
        
        % COMPUTE MACRO-SCALE STRAINS       
        % Pre-allocation to store the macro-scale Gauss Point (GP) strains
        strain = cell(1,Mesh.noElements);
        
        % Compute the macro-scale GP strains
        for idx = 1:Mesh.noElements
            
           % Set-up local node coordinates and connectivity
            connect_element = Mesh.connect(:,idx);
            x_element       = Mesh.x(:,connect_element);

            % Get the local displacement and phase-field
            switch MacroSpaceDim
                
                case 1
                    uDofMap           = connect_element;                    
                    localU            = Sol(uDofMap);
                    
                    % Loop over integration points
                    for gauss_point = gp

                        % Compute generic shape functions
                        [N, dN, ~] = ShapeFunction(x_element, gauss_point, ElementData.type);

                        % Compute shape functions for displacement and phase-field
                        [~, ~, ~, Bs] = ShapeFunctionMultiDim(N, dN, ProblemData, ElementData);

                        % Get strain, pf and its gradient at integration point
                        % (gauss_point)
                        strain{1,idx} = Bs * localU;
                    end
                    
                case 2           
                    uDofMap           = zeros(n,1);
                    uDofMap(1:2:end)  = 2*connect_element-1;
                    uDofMap(2:2:end)  = 2*connect_element;
                    localU            = Sol(uDofMap);

                    % Loop over integration points
                    for gauss_point = gp

                        % Compute generic shape functions
                        [N, dN, ~] = ShapeFunction(x_element, gauss_point, ElementData.type);

                        % Compute shape functions for displacement and phase-field
                        [~, BU, ~, ~] = ShapeFunctionMultiDim(N, dN, ProblemData, ElementData);

                        % Get strain, pf and its gradient at integration point
                        % (gauss_point)
                        strain{1,idx} = BU * localU;

                    end
                case 3
                    error('3D not yet implemented!')
                otherwise
                    error('Wrong input!')
            end                   
        end
        
        % SOLVE THE RVE PROBLEMS
        % Pre-allocation to store the macro-scale Gauss Point (GP) stresses
        stress = cell(1,Mesh.noElements);
        D      = cell(1,Mesh.noElements);
        
        parfor idx = 1:Mesh.noElements
            t = getCurrentTask();
            taskID=t.ID;
            try
                mphstart(comsolPort(taskID));
            catch
            end 
            import('com.comsol.model.*');
            import('com.comsol.model.util.*');
            
            [D{1,idx},stress{1,idx},flag{1,idx}] = RVESolve(idx,strain{1,idx},...
                                         macroSpaceDim,'Jac',RVEdata);
        end
        
        % Store the homogenised stresses
        GPdata.stress = stress;
        
        % ASSEMBLE MACRO-SCALE RESIDUAL & STIFFNESS MATRIX
        
        % Check if all RVE problems were successfully solved
        for idx = 1:Mesh.noElements
            el_flag = flag{1,idx};
            if el_flag(1) == 0
                jac_flag = 0;
                break;
            else
                jac_flag = 1;
            end
            if el_flag(2) == 0
                res_flag = 0;
                break;
            else
                res_flag = 1;
            end
        end
        
        % Assemble only if all RVE problems were solved successfully
        if res_flag == 1
            
            % Pre-allocate the macro-scale residual vector
            Res = zeros(NN,1);
            
            % Pre-allocate the macro-scale stiffness matrix
            if jac_flag == 1
                K = sparse([], [], [], NN, NN, 20*NN);
            end
        
            % Loop over elements
            for idx = 1:Mesh.noElements

                % Set-up local node coordinates and connectivity
                connect_element = Mesh.connect(:,idx);
                x_element       = Mesh.x(:,connect_element);

                switch macroSpaceDim

                    case 1 % 1D
                        uDofMap           = connect_element;

                        % Loop over integration points
                        for gauss_point = gp
                            % Compute generic shape functions
                            [N, dN, j] = ShapeFunction(x_element, gauss_point, MacroData.problem.type);

                            % Compute shape functions for displacement and phase-field
                            [~, ~, ~, Bs] = ShapeFunctionMultiDim(N, dN, MacroData.problem);

                            % Reference to real area mapping
                            dX = gauss_point(1) * j;

                            % Compute the residual (displacement)
                            Res_e = Bs' * stress{1,idx} * dX * MacroData.problem.thickness{1,idx};
                            
                            if jac_flag == 1 
                                K_e = Bs' * D{1,idx} * Bs * dX * MacroData.problem.thickness{1,idx};
                            end
                                

                        end

                    case 2 % 2D               
                        % Get the local displacement and phase-field (2-D)
                        uDofMap           = zeros(n,1);
                        uDofMap(1:2:end)  = 2*connect_element-1;
                        uDofMap(2:2:end)  = 2*connect_element;

                        % Loop over integration points
                        for gauss_point = gp
                            % Compute generic shape functions
                            [N, dN, j] = ShapeFunction(x_element, gauss_point, MacroData.problem.type);

                            % Compute shape functions for displacement and phase-field
                            [~, BU, ~, ~] = ShapeFunctionMultiDim(N, dN, MacroData.problem);

                            % Reference to real area mapping
                            dX = gauss_point(1) * j;

                            % Compute the residual (displacement)
                            Res_e = BU' * stress{1,idx} * dX * MacroData.problem.thickness{1,idx};
                            
                            if jac_flag == 1 
                                K_e = BU' * D{1,idx} * BU * dX * MacroData.problem.thickness{1,idx};
                            end

                        end

                    case 3
                        error('3D not yet implemented!')
                    otherwise
                        error('Wrong input!')
                end                    

                % Assemble the global residuals
                Res (uDofMap)  = Res (uDofMap)  + Res_e;
                
                if jac_flag == 1
                    K (uDofMap,uDofMap)  = K (uDofMap,uDofMap)  + K_e;
                end

            end

            % Check dimensions of the assembled macro-scale residual
            assert(length(Res)  == Mesh.noNodes * MacroData.problem.dof, 'Incorrect Res dimensions!');
            
            if jac_flag == 1
                % Check dimensions of the assembled macro-scale stiffness matrix
                assert(length(K(:,1))  == Mesh.noNodes * ProblemData.dof, 'Incorrect K dimensions!');
                assert(length(K(1,:))  == Mesh.noNodes * ProblemData.dof, 'Incorrect K dimensions!');
            end

            K       = sparse(K);
            Mflag(2) = res_flag;
            Mflag(1) = jac_flag;
            
        elseif res_flag == 0
            
            Res = [];
            K   = [];
            Mflag(2) = res_flag;
            Mflag(1) = 0;
            
        else            
            error('Check .\src\macro\AssembleMacroSystem.m for errors!')
            
        end     
        
    otherwise
        error('Wrong input!')
end
            
end

