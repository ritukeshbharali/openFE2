%-------------------------------------------------------------------------%
% [D,stress,flag] = RVESolvePFDAM1(idx,strain,MacroSpaceDim,action)
% solves the RVE problem for element 'idx', using the macro-scale 'strain'
% and yields the macro-scale/homogenised 'stress' and/or the homogenised
% tangent stiffness moduli 'D' via sensitivity analysis.
%
% INPUT:  idx            -> Macro-scale Element
%         strain         -> Macro-scale strain (Voigt notation)
%         MacroSpaceDim  -> Spatial Dimension of Macro-scale problem
%         action         -> 'Jac' if D-tensor computation required
%
% OUTPUT: D              -> Macro-scale homogenised tangent moduli
%         stress         -> Macro-scale stress (Voigt notation)
%         flag           -> Indicates whether RVE computation was a success
%                           (1) or failure (0). Format: [D-tensor,
%                                                        RVE Problem]
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   06.01.2021
%-------------------------------------------------------------------------%

function [D,stress,flag] = RVESolvePFDAM1(idx,strain,MacroSpaceDim,...
                                             action)

% Load the RVE corresponding to the element 'idx'
RVEname = sprintf('RVE%d.mph', idx);
folder = './output/RVEs';
model = mphload(fullfile(folder,RVEname));

% Set flag values [D-tensor, Residual] to 1
flag = [1,1];

switch MacroSpaceDim
    
    % 1D macro-scale problem
    case 1
        
        % Get strain to be imposed on the RVE
        model.param.set('eps_xx', strain(1));
        model.param.set('eps_yy', 0);
        model.param.set('eps_xy', 0);
        
        try
            % Solve the RVE problem
            model.sol('sol1').runAll;

            % Compute the homogenised dual quantities (stress)
            model.result.numerical('av1').setResult;
            tabular_data = mphtable(model,'avgtable');
            stress       = tabular_data.data(2,2);   % stress_xx for 1D
            model.result.table('avgtable').clearTableData() % Clear table            
        catch
            % Failed to solve RVE problem
            disp(['RVE ',num2str(idx),' failed to converged!']);
            stress    = [];
            flag(2)   = 0;                      
        end
        
        % Compute the homogenised D-tensor using sensitivity        
        switch action
            
            % Compute the D-tensor 
            case 'Jac'
                
                try                    
                    % Increment strain by a small amount
                    model.param.set('eps_xx', 1.001 * strain(1));
                    model.param.set('eps_yy', 0);
                    model.param.set('eps_xy', 0);
                    
                    % Solve a step for the 'small' strain increment
                    model.sol('sol2').runAll;

                    % Compute the homogenised quantities
                    model.result.numerical('av2').setResult;
                    tabular_data = mphtable(model,'avgtable2');
                    stress2 = tabular_data.data(2,2);
                    model.result.table('avgtable2').clearTableData()
                    model.sol('sol2').clearSolutionData();

                    % Compute the D-tensor (1,1)
                    D = (stress2-stress)/(1e-3*strain(1));
                catch
                    % Failed to compute the D-tensor
                    D       = [];
                    flag(1) = 0;
                end
                
            otherwise
                D       = [];
                flag(1) = 0;
                
        end
    
    % 2D macro-scale problem    
    case 2
        
        % Get strain to be imposed on the RVE
        model.param.set('eps_xx', strain(1));
        model.param.set('eps_yy', strain(2));
        model.param.set('eps_xy', strain(3));
        
        try
            % Solve the RVE problem
            model.sol('sol1').runAll;

            % Compute the homogenised dual quantities (stress)
            model.result.numerical('av1').setResult;
            tabular_data = mphtable(model,'avgtable');
            stress       = tabular_data.data(2,2:4)';   % Stress for 2D
            model.result.table('avgtable').clearTableData() % Clear table            
        catch
            % Failed to solve RVE problem
            disp(['RVE ',num2str(idx),' failed to converged!']);
            stress    = [];
            flag(2)   = 0;                      
        end
        
        % Compute the homogenised D-tensor using sensitivity        
        switch action
            
            % Compute the D-tensor 
            case 'Jac'
                
                try
                    D = zeros(3);
                    
                    % Increment strain_xx by a small amount
                    model.param.set('eps_xx', 1.001 * strain(1));
                    model.param.set('eps_yy', strain(2));
                    model.param.set('eps_xy', strain(3));
                    
                    % Solve a step for the 'small' strain increment
                    model.sol('sol2').runAll;

                    % Compute the homogenised quantities
                    model.result.numerical('av2').setResult;
                    tabular_data = mphtable(model,'avgtable2');
                    stress2 = tabular_data.data(2,2:4)';
                    model.result.table('avgtable2').clearTableData()
                    model.sol('sol2').clearSolutionData();

                    % Compute the D-tensor (First column )
                    D(:,1) = (stress2-stress)/(1e-3*strain(1));
                    
                    % Increment strain_yy by a small amount
                    model.param.set('eps_xx', strain(1));
                    model.param.set('eps_yy', 1.001 * strain(2));
                    model.param.set('eps_xy', strain(3));
                    
                    % Solve a step for the 'small' strain increment
                    model.sol('sol2').runAll;

                    % Compute the homogenised quantities
                    model.result.numerical('av2').setResult;
                    tabular_data = mphtable(model,'avgtable2');
                    stress2 = tabular_data.data(2,2:4)';
                    model.result.table('avgtable2').clearTableData()
                    model.sol('sol2').clearSolutionData();

                    % Compute the D-tensor (Second column )
                    D(:,2) = (stress2-stress)/(1e-3*strain(1));
                    
                    % Increment strain_xy by a small amount
                    model.param.set('eps_xx', strain(1));
                    model.param.set('eps_yy', strain(2));
                    model.param.set('eps_xy', 1.001 * strain(3));
                    
                    % Solve a step for the 'small' strain increment
                    model.sol('sol2').runAll;

                    % Compute the homogenised quantities
                    model.result.numerical('av2').setResult;
                    tabular_data = mphtable(model,'avgtable2');
                    stress2 = tabular_data.data(2,2:4)';
                    model.result.table('avgtable2').clearTableData()
                    model.sol('sol2').clearSolutionData();

                    % Compute the D-tensor (Third column )
                    D(:,3) = (stress2-stress)/(1e-3*strain(1));                       
                catch
                    % Failed to compute the D-tensor
                    D       = [];
                    flag(1) = 0;
                end
                
            otherwise
                D       = [];
                flag(1) = 0;
                
        end
    
    % 3D macro-scale problem   
    case 3
        error('3D RVE not implemented! Choose 1D or 2D!')
    otherwise
        error('Wrong input!')
end

% Save the model and exit safely
mphsave(model,fullfile(folder,RVEname));
end

