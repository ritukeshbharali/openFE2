%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% openFE2 is a software framework capable to carrying out multi-scale FE2
% simulations. The MACRO-SCALE problem is solved in MATLAB while the RVE
% problems are solved in COMSOL Multiphysics, in parallel. The parallel
% setup is achieved using the MATLAB Parallel Computing Toolbox.
% 
% Author: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%
% Features/Limitations: 
%          1. Bar, triangle and tetrahedon elements for 1D/2D/3D problems
%          2. Single Gauss point (makes it easier to link RVEs)
%          3. Elasticity, Phase-field Damage RVE models
%          4. Flexibility choice of RVEs solved in parallel and number of
%             cores alloted for each RVE
%          5. Restart simulation from any saved state
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clear workspace
clearvars; close all; clc
format long;

% Start timer
tic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SOFTWARE INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ') 
disp('%----------------------------------------------------------------------%')
disp('       ******* openFE2: Multi-scale FE2 analysis *******');
disp(' ')
disp('               LICENSE: openFE2 is a free software. You can')
disp('                        distribute and/or modify it as per the')
disp('                        MIT License.')
disp(' ')
disp('               AUTHOR:  Ritukesh Bharali');
disp('                        Material and Computational Mechanics Division,');
disp('                        Dept. of Industrial and Material Science,');
disp('                        Chalmers University of Technology,');
disp('                        Gothenburg, Sweden');
disp('                        mail: ritukesh.bharali@chalmers.se');
disp('%----------------------------------------------------------------------%')
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD PATH TO DIRECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('    - Setting up directories')
dir_input  = '.\input';
dir_src    = '.\src';
dir_output = '.\output';
if ~exist(dir_output, 'dir')
    mkdir(dir_output)
end
% dir_mlilib  = 'C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli';
dir_mlilib  = 'C:\Progs\COMSOL55\mli';

addpath(genpath(dir_input));
addpath(genpath(dir_src));
addpath(genpath(dir_output));
addpath(genpath(dir_mlilib));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ INPUT AND SETUP RUNTIME PARAMETERS, MACRO MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('    - Reading input data')
input_parallel
input_macro
input_rve

disp('    - Setting up macro-scale mesh')
MacroMesh = read_mesh(dir_input,inp_macro.problem);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP COMSOL SERVERS FOR RVE PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('    - Setting up COMSOL servers')
comsolPort = linspace(2036,2036+inp_parallel.nparRVE-1,inp_parallel.nparRVE);

% Start COMSOL mphservers
for num = 1:inp_parallel.nparRVE
    COMSOLPort = comsolPort(num);
    % Add path to comsolmphserver.exe below within " "
    % system( ['"C:\Program Files\COMSOL\COMSOL55\Multiphysics\bin\win64\comsolmphserver.exe" -np ' num2str(inp_parallel.ncoreRVE) ' -port ' num2str(COMSOLPort) ' &'] );
    system( ['"C:\Progs\COMSOL55\bin\win64\comsolmphserver.exe" -np ' num2str(inp_parallel.ncoreRVE) ' -port ' num2str(COMSOLPort) ' &'] );
    pause(.1)
end
clear num COMSOLPort;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALLOCATE RVEs TO EACH MACRO GAUSS POINT AND COMPUTE INITIAL STIFFNESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('    - Allocating RVEs & computing initial stiffness matrix')
[K,D,RVEdata] = AssembleInitMacroSystem(MacroMesh,...
                             inp_macro.problem,inp_rve.problem,comsolPort);

                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP A GLOBAL DATABASE AND STORE HOMOGENISED TANGENT MODULI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('    - Setting up a global database')
globDat    = setupGlobalDatabase(K,MacroMesh,inp_macro,RVEdata);
globDat.D  = D;
clear D;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME-STEPPING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise vectors to store solution and residuals
NN   = MacroMesh.noNodes * inp_macro.problem.dof;
Sol  = zeros(NN,1);
Res  = zeros(NN,1);

% Construct a constraint matrix [Shephard (1984)]
CMat = speye(NN);
CMat(:,inp_macro.constraints(:,1)) = [];

% Begin time stepping here
step     = 1;
time     = 0.0;
timestep = inp_macro.solver.initstepsize;

% Print information to the command window
fprintf('\n')
disp('-------------------------------------------------------')
disp(['Step ', num2str(step), ': Time = ', num2str(time)])
disp('-------------------------------------------------------')

disp(['Reaction = ',num2str(globDat.Force(inp_macro.postproc.loadDof,step)), ...
      '[N], Displacement = ',num2str(globDat.Disp(inp_macro.postproc.dispDof,step)), ' [m]'])

while true
    
    % Increment step
    step = step + 1;
    
    % Step-size cut
    if ismember(step, inp_macro.solver.cutsteps) == 1
        timestep = timestep/inp_macro.solver.cutstepsize;
    end
    
    % Time increment
    time = time + timestep;
    
    % Print information to the command window
    fprintf('\n')
    disp('-------------------------------------------------------')
    disp(['Step ', num2str(step), ': Time = ', num2str(time)])
    disp('-------------------------------------------------------')
    
    % Set iteration to zero
    iter = 0;
    
    % Apply BC
    Sol = ApplyDBC(Sol,time,inp_macro.constraints);
    
    % Store iteration 0 solution
    Sol0 = Sol;
    
    % Begin Newton iterations
    while true
        
        % Iteration counter increment
        iter = iter + 1;
        
        % Store a copy of old solutions, residuals, jacobians
        % Previous iteration values       
        oldSol   = Sol;
        oldRes   = Res;
        oldK     = K;
        
        % Modify history variable (Comment for elasticity!)
         RVESetHistory(MacroMesh,inp_rve.problem.type,...
                                       globDat.RVEoldStepSol,comsolPort);
        
        % Get macro-scale residual and/or stiffness matrix        
        [K,Res,GPData,Mflag] = AssembleMacroSystem(MacroMesh,Sol,...
                                         inp_macro,RVEdata,...
                                         comsolPort,step,iter);
        
        % If stiffness matrix is not returned from AssembleMacroSystem, use
        % the old stiffness matrix       
        if Mflag(1) == 0
            K = oldK;
        end
        
        % If residual is not returned from AssembleMacroSystem or the
        % maximum number of iterations are reached, revert to the previous
        % converged step and reduce the timestep size.
        if Mflag(2) == 0 || iter == inp_macro.solver.maxiter
            
            % Display message to command window
			disp('Macro-problem failed to converge!');
            
            % Revert macro-scale variables
			Sol   = globDat.oldStepSol;
			Res   = globDat.oldStepRes;
			K     = globDat.oldStepK;
			time  = time - timestep;
			step  = step - 1;
			
			% Reduce timestep size
			timestep = timestep/inp_macro.solver.cutstepsize;
            
            % Revert all RVEs to old 'converged' step
            RVERevert(MacroMesh,inp_rve.problem.type,...
                                      globDat.RVEoldStepSol,comsolPort);
            
            % Display messages to the command window
			disp(['Reverted to Step',num2str(step)]);
            disp(['Step-size set to ',num2str(timestep)]);
            
            break;
        end
        
        % Solve the condensed linear problem
        [SolUpdate,LinSolverData] = LinearSolver(CMat'*K*CMat,CMat'*Res,...
                                                inp_macro.solver);
               
        % Update solution
        SolUpdate = CMat*SolUpdate;
        Sol       = Sol - SolUpdate;
        
        % Compute error in solution field
        err = norm(Sol-oldSol)/norm(Sol0);
        
        % Display errors in the command window
		disp(['Iter = ',num2str(iter),', Error in Sol = ',num2str(err)])
        
        % Check for convergence        
        if err < inp_macro.solver.nonlinsolvertol
            
            % Display message to command window
            disp(['Macro-problem converged in ', ...
                                      num2str(iter),' iteration(s).'])
                                  
            % Modify history variable (Comment for elasticity!)
             RVESetHistory(MacroMesh,inp_rve.problem.type,...
                                       globDat.RVEoldStepSol,comsolPort);                                  
                                  
            [~,Res,GPData,Mflag] = AssembleMacroSystem(MacroMesh,Sol,...
                                         inp_macro,RVEdata,...
                                         comsolPort,'NoJac');
                                     
            % If residual is not returned from AssembleMacroSystem,
            % revert to previous converged step and reduce the timestep size.
            if Mflag(2) == 0

                % Display 'failure to converge' message to command window
                disp('Internal forces computation failed!');

                % Revert macro-scale variables
                Sol   = globDat.oldStepSol;
                Res   = globDat.oldStepRes;
                K     = globDat.oldStepK;
                time  = time - timestep;
                step  = step - 1;

                % Reduce timestep size
                timestep = timestep/inp_macro.solver.cutstepsize;

                % Revert all RVEs to old 'converged' step
                RVERevert(MacroMesh,inp_rve.problem.type,...
                                          globDat.RVEoldStepSol,comsolPort);

                % Display messages to the command window
                disp(['Reverted to Step',num2str(step)]);
                disp(['Step-size set to ',num2str(timestep)]);

                break;
            end            
                                     
            % Store nodal forces and solution into global database
			globDat.Force(:,step) = Res;
			globDat.Disp(:,step)  = Sol;
            globDat.stress        = GPData.stress;
            
            % Update oldoldStep data
			globDat.oldoldStepSol  = globDat.oldStepSol;
			globDat.oldoldStepRes  = globDat.oldStepRes;
			globDat.oldoldStepK    = globDat.oldStepK;
            
            % Update oldStep data
			globDat.oldStepSol  = Sol;
			globDat.oldStepRes  = Res;
			globDat.oldStepK    = K;
            
            disp(['Reaction = ',num2str(globDat.Force(inp_macro.postproc.loadDof,step)), ...
                      '[N], Displacement = ',num2str(globDat.Disp(inp_macro.postproc.dispDof,step)), ' [m]'])
                         
            % Print out Macro-GPData and RVE images
            [globDat.RVEoldStepSol,~] = RVEPostProcess(MacroMesh,...
                    step,inp_rve.problem,comsolPort);
                
            % Create restart files if required
            if ismember(step,inp_macro.restart.steps) == 1
                CreateRestartFiles(Sol,Res,K,MacroMesh,RVEdata,step,time);
            end
            
            break;
            
        end   
       
    % End Newton iterations    
    end
    
    % Stop timestepping if certain criteria is met
    if step     == inp_macro.solver.nsteps || ...
       timestep == inp_macro.solver.finalstepsize
       disp(' ')
       disp('    - Analysis termination criterion reached!')
       break;
    end

% End time-stepping    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START ADDITIONAL POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a folder to store post-processing stuff if it does not exist yet
dir_postproc = './output/images';
if ~exist(dir_postproc, 'dir')
    mkdir(dir_postproc)
end

% Resize globDat Force and Solution
globDat.Force(:,step+1:end) = [];
globDat.Disp(:,step+1:end) = [];

figure
plot(globDat.Disp(inp_macro.postproc.dispDof,:),...
    -globDat.Force(inp_macro.postproc.loadDof,:),'-k','LineWidth',1.5)
hold on;
xlabel('Macro-displacement [m]')
ylabel('Macro-reaction [N]')
set(gca,'FontSize', 15);
saveas(gcf,fullfile(dir_postproc,'lodi.fig'))
saveas(gcf,fullfile(dir_postproc,'lodi.png'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('*** This is the end, my friend, the very end! ***')
disp(['*** TIME ELAPSED = ',num2str(toc)])
disp(' ')

% Kill COMSOLmphserver processes (MAY NOT BE SAFE!)
!taskkill -im cmd.exe 