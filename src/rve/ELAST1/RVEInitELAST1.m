%-------------------------------------------------------------------------%
% [D,RVEdofs] = RVEInitELAST1(idx,MacroSpaceDim,RVEdata) initialises a
% linear elastic RVE with an inclusion and computes the homogenised tangent 
% stiffness moduli 'D'.
%
% INPUT:  idx            -> Macro-scale Element
%         MacroSpaceDim  -> Spatial Dimension of Macro-scale problem
%         RVEdata        -> MATLAB struct, must contain the spatial
%                           dimension of the RVE problem 'SpaceDim'
%
% OUTPUT: D              -> Macro-scale homogenised tangent moduli
%         RVEdofs        -> RVE number of degrees of freedom
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   06.01.2021
%
%-------------------------------------------------------------------------%

function [D,RVEdofs] = RVEInitELAST1(idx,MacroSpaceDim,RVEdata)

% Import COMSOL libraries if not done!
import com.comsol.model.*
import com.comsol.model.util.*

% Give the RVE file a name
RVEname = sprintf('RVE%d.mph', idx);

% Create a folder to store the RVEs if it does not exist
folder = './output/RVEs';
if ~exist(folder, 'dir')
    mkdir(folder)
end

% Delete if exisitng model with same filename is found!
if isfile(fullfile(folder,RVEname))
    delete(fullfile(folder,RVEname));
    disp(['     -> Deleted ', RVEname])
else
    disp(['     -> Initialising RVE ',num2str(idx)])
end

switch RVEdata.SpaceDim
    
    % 1D RVE
    case 1
        error('1D RVE not implemented!')
    
    % 2D RVE
    case 2
        
        %    -------
        %    |  o  |
        %    -------
        %
        %    RVE with stiff inclusion
        
        
        % Create a 2D RVE model
        model = ModelUtil.create('Model');
        
        % Uncomment if you wish to observe model progress in separate
        % windows. Note that one window will pop up for each process!
        % ModelUtil.showProgress(true);
        
        % Set runtime parameters (Modify as required!)
        model.param.set('E_incl', '210[GPa]');
        model.param.set('E_mat', '100[GPa]');
        model.param.set('nu_incl', '0.4');
        model.param.set('nu_mat', '0.25');
        model.param.set('hmax', '3e-2[mm]');
        model.param.set('eps_xx', '1e-10');
        model.param.set('eps_yy', '0');
        model.param.set('eps_xy', '0');
        
        % Create component, geometry and mesh
        model.component.create('comp1', true);
        model.component('comp1').geom.create('geom1', 2);
        model.component('comp1').geom('geom1').lengthUnit('mm');
        model.component('comp1').geom('geom1').useConstrDim(false);
        model.component('comp1').geom('geom1').create('sq1', 'Square');
        model.component('comp1').geom('geom1').feature('sq1').set('base', 'center');
        model.component('comp1').geom('geom1').create('c1', 'Circle');
        model.component('comp1').geom('geom1').feature('c1').set('r', 0.25);
        model.component('comp1').geom('geom1').run;
        model.component('comp1').sorder('linear');
        model.component('comp1').curvedInterior(false);
        model.component('comp1').mesh.create('mesh1');
        model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
        model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
        model.component('comp1').mesh('mesh1').feature('size').set('hmax', 'hmax');
        model.component('comp1').mesh('mesh1').feature('size').set('hmin', 'hmax');
        model.component('comp1').mesh('mesh1').run;
        
        % Some additional definitions for periodicity conditions
        model.component('comp1').cpl.create('linext1', 'LinearExtrusion');
        model.component('comp1').cpl.create('linext2', 'LinearExtrusion');
        model.component('comp1').cpl('linext1').selection.geom('geom1', 1);
        model.component('comp1').cpl('linext1').selection.set([1]);
        model.component('comp1').cpl('linext2').selection.geom('geom1', 1);
        model.component('comp1').cpl('linext2').selection.set([2]);
        model.component('comp1').cpl('linext1').label('mapping x');
        model.component('comp1').cpl('linext1').set('opname', 'linextx');
        model.component('comp1').cpl('linext1').set('method', 'closest');
        model.component('comp1').cpl('linext1').selection('srcvertex1').set([1]);
        model.component('comp1').cpl('linext1').selection('srcvertex2').set([2]);
        model.component('comp1').cpl('linext1').selection('dstvertex1').set([7]);
        model.component('comp1').cpl('linext1').selection('dstvertex2').set([8]);
        model.component('comp1').cpl('linext2').label('mapping y');
        model.component('comp1').cpl('linext2').set('opname', 'linexty');
        model.component('comp1').cpl('linext2').set('method', 'closest');
        model.component('comp1').cpl('linext2').selection('srcvertex1').set([1]);
        model.component('comp1').cpl('linext2').selection('srcvertex2').set([7]);
        model.component('comp1').cpl('linext2').selection('dstvertex1').set([2]);
        model.component('comp1').cpl('linext2').selection('dstvertex2').set([8]);

        % Periodicity in x-direction
        model.component('comp1').variable.create('var1');
        model.component('comp1').variable('var1').set('dx', 'x-linextx(x)');
        model.component('comp1').variable('var1').set('dy', 'y-linextx(y)');
        model.component('comp1').variable('var1').set('du', 'u-linextx(u)');
        model.component('comp1').variable('var1').set('dv', 'v-linextx(v)');
        model.component('comp1').variable('var1').selection.geom('geom1', 1);
        model.component('comp1').variable('var1').selection.set([4]);
        model.component('comp1').variable('var1').label('on positive x');

        % Periodicity in y-direction
        model.component('comp1').variable.create('var2');
        model.component('comp1').variable('var2').set('dx', 'x-linexty(x)');
        model.component('comp1').variable('var2').set('dy', 'y-linexty(y)');
        model.component('comp1').variable('var2').set('du', 'u-linexty(u)');
        model.component('comp1').variable('var2').set('dv', 'v-linexty(v)');
        model.component('comp1').variable('var2').selection.geom('geom1', 1);
        model.component('comp1').variable('var2').selection.set([3]);
        model.component('comp1').variable('var2').label('on positive y');
        
        % Solid Mechanics module (Momentum balance equation)
        model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');
        model.component('comp1').physics('solid').create('lemm2', 'LinearElasticModel', 2);
        model.component('comp1').physics('solid').feature('lemm2').selection.set([2]);
        model.component('comp1').physics('solid').create('fix1', 'Fixed', 0);
        model.component('comp1').physics('solid').feature('fix1').selection.set([1]);
        model.component('comp1').physics('solid').create('constr1', 'PointwiseConstraint', 1);
        model.component('comp1').physics('solid').feature('constr1').selection.set([3 4]);
        model.component('comp1').physics('solid').create('constr2', 'PointwiseConstraint', 1);
        model.component('comp1').physics('solid').feature('constr2').selection.set([3 4]);
        model.component('comp1').physics('solid').prop('ShapeProperty').set('order_displacement', 1);
        model.component('comp1').physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 'Quasistatic');
        model.component('comp1').physics('solid').feature('lemm1').set('E_mat', 'userdef');
        model.component('comp1').physics('solid').feature('lemm1').set('E', 'E_mat');
        model.component('comp1').physics('solid').feature('lemm1').set('nu_mat', 'userdef');
        model.component('comp1').physics('solid').feature('lemm1').set('nu', 'nu_mat');
        model.component('comp1').physics('solid').feature('lemm2').set('E_mat', 'userdef');
        model.component('comp1').physics('solid').feature('lemm2').set('E', 'E_incl');
        model.component('comp1').physics('solid').feature('lemm2').set('nu_mat', 'userdef');
        model.component('comp1').physics('solid').feature('lemm2').set('nu', 'nu_incl');
        model.component('comp1').physics('solid').feature('constr1').set('constraintExpression', '-du+eps_xx*dx+eps_xy*dy');
        model.component('comp1').physics('solid').feature('constr1').set('shapeFunctionType', 'shgp');
        model.component('comp1').physics('solid').feature('constr1').set('order', 4);
        model.component('comp1').physics('solid').feature('constr1').set('frame', 'spatial');
        model.component('comp1').physics('solid').feature('constr2').set('constraintExpression', '-dv+eps_xy*dx+eps_yy*dy');
        model.component('comp1').physics('solid').feature('constr2').set('shapeFunctionType', 'shgp');
        model.component('comp1').physics('solid').feature('constr2').set('order', 4);
        model.component('comp1').physics('solid').feature('constr2').set('frame', 'spatial');
        
        % Define Solver
        model.study.create('std2');
        model.study('std2').create('time', 'Transient');
        model.sol.create('sol1');
        model.sol('sol1').study('std2');
        model.sol('sol1').attach('std2');
        model.sol('sol1').create('st1', 'StudyStep');
        model.sol('sol1').create('v1', 'Variables');
        model.sol('sol1').create('t1', 'Time');
        model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
        model.sol('sol1').feature('t1').feature.remove('fcDef');
        model.study('std2').feature('time').set('tlist', 'range(0,0.5,1)');
        model.study('std2').feature('time').set('usertol', true);
        model.study('std2').feature('time').set('rtol', '1e-3');
        model.study('std2').feature('time').set('useinitsol', true);
        model.study('std2').feature('time').set('usesol', true);
        model.sol('sol1').attach('std2');
        model.sol('sol1').feature('v1').set('clist', {'range(0,0.5,1)' '0.001[s]'});
        model.sol('sol1').feature('v1').feature('comp1_u').set('scalemethod', 'manual');
        model.sol('sol1').feature('v1').feature('comp1_u').set('scaleval', '1e-2*0.0014142135623730952');
        model.sol('sol1').feature('t1').set('tlist', 'range(0,0.5,1)');
        model.sol('sol1').feature('t1').set('rtol', '1e-3');
        model.sol('sol1').feature('t1').set('timemethod', 'genalpha');
        model.sol('sol1').feature('t1').set('tstepsgenalpha', 'strict');
        model.sol('sol1').feature('t1').set('initialstepgenalpha', 1);
        model.sol('sol1').feature('t1').set('initialstepgenalphaactive', true);
        model.sol('sol1').feature('t1').set('maxstepconstraintgenalpha', 'const');
        model.sol('sol1').feature('t1').set('maxstepgenalpha', 1);
        model.sol('sol1').feature('t1').set('rhoinf', 0.85);
        model.sol('sol1').feature('t1').set('rescaleafterinitbw', true);
        model.sol('sol1').feature('t1').feature('dDef').set('ooc', false);
        model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
        model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 25);
        model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'onevery');
        
        % Create post-processing step to compute homogenised dual
        % quantities
        model.result.table.create('avgtable', 'Table');
        model.result.numerical.create('av1', 'AvSurface');
        model.result.numerical('av1').set('data', 'dset1');
        model.result.numerical('av1').selection.all;
        model.result.numerical('av1').set('probetag', 'none');
        model.result.numerical('av1').set('table', 'avgtable');
        model.result.numerical('av1').set('expr', {'solid.sx' 'solid.sy' 'solid.sxy'});
        model.result.numerical('av1').set('unit', {'N/m^2' 'N/m^2' 'N/m^2'});
        model.result.numerical('av1').set('descr', {'Stress tensor, x component' 'Stress tensor, y component' 'Stress tensor, xy component'});
        model.result.numerical('av1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
        
        % Create a surface plot for visualisation
        model.result.create('pg1', 'PlotGroup2D');        
        model.result('pg1').set('data', 'dset1');    
        model.result('pg1').create('surf1', 'Surface');
        
        % Compute the D-tensor based on spatial dimension of the
        % macro-scale problem
        switch MacroSpaceDim
            
            % 1D macro-scale problem
            case 1
                
                % Set a small strain in x-direction
                model.param.set('eps_xx', '1e-10');
                model.param.set('eps_yy', '0');
                model.param.set('eps_xy', '0');
                
                % Solve the RVE problem
                model.sol('sol1').runAll;
                
                % Compute homogenised stress and sensitivity
                model.result.numerical('av1').setResult;
                tabular_data = mphtable(model,'avgtable');
                stress       = tabular_data.data(2,2);
                D = stress/1e-10;
                
                % Clear table
                model.result.table('avgtable').clearTableData();
                
                % Extract some additional data from COMSOL
                solinfo         = mphsolinfo(model,'soltag','sol1');
                RVEdofs         = solinfo.size;
                
                % Change solver settings for subsequent steps
                model.sol('sol1').feature('t1').set('initialstepgenalphaactive', true);
                model.sol('sol1').feature('t1').set('initialstepgenalpha', 1);
                model.sol('sol1').feature('t1').set('tstepsgenalpha', 'strict');
                model.study('std2').feature('time').set('useinitsol', true);
                model.study('std2').feature('time').set('initmethod', 'sol');
                model.study('std2').feature('time').set('initstudy', 'std2');
                model.study('std2').feature('time').set('solnum', 2);
                model.sol('sol1').feature('v1').set('initmethod', 'sol');
                model.sol('sol1').feature('v1').set('initsol', 'sol1');
                model.sol('sol1').feature('v1').set('solnum', 2);
                
                % Create a second study (Required to compute D-tensor via
                % sensitivity analysis)
                model.study.create('std3');
                model.study('std3').create('time', 'Transient');
                model.sol.create('sol2');
                model.sol('sol2').study('std3');
                model.sol('sol2').attach('std3');
                model.sol('sol2').create('st1', 'StudyStep');
                model.sol('sol2').create('v1', 'Variables');
                model.sol('sol2').create('t1', 'Time');
                model.sol('sol2').feature('t1').create('fc1', 'FullyCoupled');
                model.sol('sol2').feature('t1').create('st1', 'StopCondition');
                model.sol('sol2').feature('t1').feature.remove('fcDef');
                model.study('std3').feature('time').set('tlist', 'range(0,0.5,1)');
                model.study('std3').feature('time').set('usertol', true);
                model.study('std3').feature('time').set('rtol', '1e-3');
                model.study('std3').feature('time').set('useinitsol', true);
                model.study('std3').feature('time').set('initmethod', 'sol');
                model.study('std3').feature('time').set('initstudy', 'std2');
                model.study('std3').feature('time').set('solnum', 2);
                model.sol('sol2').feature('v1').set('initmethod', 'sol');
                model.sol('sol2').feature('v1').set('initsol', 'sol1');
                model.sol('sol2').feature('v1').set('solnum', 2);
                model.sol('sol2').attach('std3');
                model.sol('sol1').feature('v1').set('clist', {'range(0,0.5,1)' '0.001[s]'});
                model.sol('sol1').feature('v1').feature('comp1_u').set('scalemethod', 'manual');
                model.sol('sol1').feature('v1').feature('comp1_u').set('scaleval', '1e-2*0.0014142135623730952');
                model.sol('sol1').feature('t1').set('tlist', 'range(0,0.5,1)');
                model.sol('sol1').feature('t1').set('rtol', '1e-3');
                model.sol('sol1').feature('t1').set('timemethod', 'genalpha');
                model.sol('sol1').feature('t1').set('tstepsgenalpha', 'strict');
                model.sol('sol1').feature('t1').set('initialstepgenalpha', 1);
                model.sol('sol1').feature('t1').set('initialstepgenalphaactive', true);
                model.sol('sol1').feature('t1').set('maxstepconstraintgenalpha', 'const');
                model.sol('sol1').feature('t1').set('maxstepgenalpha', 1);
                model.sol('sol1').feature('t1').set('rhoinf', 0.85);
                model.sol('sol1').feature('t1').set('rescaleafterinitbw', true);
                model.sol('sol1').feature('t1').feature('dDef').set('ooc', false);
                model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
                model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 25);
                model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'onevery');
                model.sol('sol2').feature('t1').feature('st1').set('stopcondterminateon', {'true'});
                model.sol('sol2').feature('t1').feature('st1').set('stopcondActive', {'on'});
                model.sol('sol2').feature('t1').feature('st1').set('stopconddesc', {'Stop expression 1'});
                model.sol('sol2').feature('t1').feature('st1').set('stopcondarr', {'timestep<1e-7[s]'});


                % Create a second table for homogenised quantities (to compute
                % sensitivity)
                model.result.table.create('avgtable2', 'Table');
                model.result.numerical.create('av2', 'AvSurface');
                model.result.numerical('av2').set('data', 'dset2');
                model.result.numerical('av2').selection.all;
                model.result.numerical('av2').set('probetag', 'none');
                model.result.numerical('av2').set('table', 'avgtable2');
                model.result.numerical('av2').set('expr', {'solid.sx' 'solid.sy' 'solid.sxy'});
                model.result.numerical('av2').set('unit', {'N/m^2' 'N/m^2' 'N/m^2'});
                model.result.numerical('av2').set('descr', {'Stress tensor, x component' 'Stress tensor, y component' 'Stress tensor, xy component'});
                model.result.numerical('av2').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});

                % Save the RVE model to folder 'RVE' and exit safely
                mphsave(model,fullfile(folder,RVEname));    
            
            % 2D macro-scale problem    
            case 2
                
                % Pre-allocate D-tensor
                D = zeros(3);

                % Set a small strain in x-direction 
                model.param.set('eps_xx', '1e-10');
                model.param.set('eps_yy', '0');
                model.param.set('eps_xy', '0');
                
                % Solve the RVE problem
                model.sol('sol1').runAll;
                
                % Compute homogenised stresses and sensitivity
                model.result.numerical('av1').setResult;
                tabular_data = mphtable(model,'avgtable');
                stress       = tabular_data.data(2,2:4);
                D(:,1) = stress'/1e-10;
                
                % Clear table
                model.result.table('avgtable').clearTableData();
                
                % Clear solution
                model.sol('sol1').clearSolutionData();

                % Set a small strain in y-direction
                model.param.set('eps_xx', '0');
                model.param.set('eps_yy', '1e-10');
                model.param.set('eps_xy', '0');
                
                % Solve the RVE problem
                model.sol('sol1').runAll;
                
                % Compute homogenised stresses and sensitivity
                model.result.numerical('av1').setResult;
                tabular_data = mphtable(model,'avgtable');
                stress = tabular_data.data(2,2:4);
                D(:,2) = stress'/1e-10;
                
                % Clear table
                model.result.table('avgtable').clearTableData();
                
                % Clear solution
                model.sol('sol1').clearSolutionData();

                % Set a small strain in the xy-direction
                model.param.set('eps_xx', '0');
                model.param.set('eps_yy', '0');
                model.param.set('eps_xy', '1e-10');
                
                % Solve the RVE problem
                model.sol('sol1').runAll;
                
                % Compute homogenised stresses and sensitivity
                model.result.numerical('av1').setResult;
                tabular_data = mphtable(model,'avgtable');
                stress = tabular_data.data(2,2:4);
                D(:,3) = stress'/1e-10;
                
                % Clear table
                model.result.table('avgtable').clearTableData();  
                
                % Extract some additional data from COMSOL
                solinfo         = mphsolinfo(model,'soltag','sol1');
                RVEdofs         = solinfo.size;
                
                % Change solver settings for subsequent steps
                model.sol('sol1').feature('t1').set('initialstepgenalphaactive', true);
                model.sol('sol1').feature('t1').set('initialstepgenalpha', 1);
                model.sol('sol1').feature('t1').set('tstepsgenalpha', 'strict');
                model.study('std2').feature('time').set('useinitsol', true);
                model.study('std2').feature('time').set('initmethod', 'sol');
                model.study('std2').feature('time').set('initstudy', 'std2');
                model.study('std2').feature('time').set('solnum', 2);
                model.sol('sol1').feature('v1').set('initmethod', 'sol');
                model.sol('sol1').feature('v1').set('initsol', 'sol1');
                model.sol('sol1').feature('v1').set('solnum', 2);
                
                % Create a second study (Required to compute D-tensor via
                % sensitivity analysis)
                model.study.create('std3');
                model.study('std3').create('time', 'Transient');
                model.sol.create('sol2');
                model.sol('sol2').study('std3');
                model.sol('sol2').attach('std3');
                model.sol('sol2').create('st1', 'StudyStep');
                model.sol('sol2').create('v1', 'Variables');
                model.sol('sol2').create('t1', 'Time');
                model.sol('sol2').feature('t1').create('se1', 'Segregated');
                model.sol('sol2').feature('t1').create('ps2', 'PreviousSolution');
                model.sol('sol2').feature('t1').create('st1', 'StopCondition');
                model.sol('sol2').feature('t1').feature('se1').create('ss1', 'SegregatedStep');
                model.sol('sol2').feature('t1').feature('se1').create('ss2', 'SegregatedStep');
                model.sol('sol2').feature('t1').feature.remove('fcDef');
                model.study('std3').feature('time').set('tlist', 'range(0,0.5,1)');
                model.study('std3').feature('time').set('usertol', true);
                model.study('std3').feature('time').set('rtol', '1e-3');
                model.study('std3').feature('time').set('useinitsol', true);
                model.study('std3').feature('time').set('initmethod', 'sol');
                model.study('std3').feature('time').set('initstudy', 'std2');
                model.study('std3').feature('time').set('solnum', 2);
                model.sol('sol2').feature('v1').set('initmethod', 'sol');
                model.sol('sol2').feature('v1').set('initsol', 'sol1');
                model.sol('sol2').feature('v1').set('solnum', 2);
                model.sol('sol2').feature('v1').set('clist', {'range(0,0.5,1)' '0.001[s]'});
                model.sol('sol2').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
                model.sol('sol2').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.0014142135623730952');
                model.sol('sol2').feature('t1').set('tlist', 'range(0,0.5,1)');
                model.sol('sol2').feature('t1').set('rtol', '1e-3');
                model.sol('sol2').feature('t1').set('timemethod', 'genalpha');
                model.sol('sol1').feature('t1').set('rhoinf', 0);
                model.sol('sol2').feature('t1').set('estrat', 'exclude');
                model.sol('sol2').feature('t1').set('maxstepconstraintgenalpha', 'const');
                model.sol('sol2').feature('t1').set('maxstepgenalpha', 1);
                model.sol('sol2').feature('t1').set('rescaleafterinitbw', true);
                model.sol('sol2').feature('t1').feature('dDef').set('ooc', false);
                model.sol('sol2').feature('t1').feature('aDef').set('cachepattern', true);
                %model.sol('sol2').feature('t1').feature('se1').set('maxsegiter', 100);
                %model.sol('sol2').feature('t1').feature('se1').set('segstabacc', 'segaacc');
                %model.sol('sol2').feature('t1').feature('se1').set('segaaccdim', 50);
                model.sol('sol2').feature('t1').feature('se1').set('segterm', 'iter');	
                model.sol('sol2').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2'});
                model.sol('sol2').feature('t1').feature('se1').feature('ss1').set('segvar', {'comp1_H'});
                model.sol('sol2').feature('t1').feature('se1').feature('ss2').set('segvar', {'comp1_u'});
                model.sol('sol2').feature('t1').feature('ps2').set('prevcomp', {'comp1_H'});
                model.sol('sol2').feature('t1').feature('st1').set('stopcondterminateon', {'true'});
                model.sol('sol2').feature('t1').feature('st1').set('stopcondActive', {'on'});
                model.sol('sol2').feature('t1').feature('st1').set('stopconddesc', {'Stop expression 1'});
                model.sol('sol2').feature('t1').feature('st1').set('stopcondarr', {'timestep<1e-7[s]'});


                % Create a second table for homogenised quantities (to compute
                % sensitivity)
                model.result.table.create('avgtable2', 'Table');
                model.result.numerical.create('av2', 'AvSurface');
                model.result.numerical('av2').set('data', 'dset2');
                model.result.numerical('av2').selection.all;
                model.result.numerical('av2').set('probetag', 'none');
                model.result.numerical('av2').set('table', 'avgtable2');
                model.result.numerical('av2').set('expr', {'solid.sx' 'solid.sy' 'solid.sxy'});
                model.result.numerical('av2').set('unit', {'N/m^2' 'N/m^2' 'N/m^2' 'N/m^2' 'N/m' 'N/m'});
                model.result.numerical('av2').set('descr', {'Stress tensor, x component' 'Stress tensor, y component' 'Stress tensor, xy component'});
                model.result.numerical('av2').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});

                % Save the RVE model to folder 'RVE' and exit safely
                mphsave(model,fullfile(folder,RVEname));                
            
            % 3D macro-scale problem    
            case 3
                error('3D RVE not implemented! Choose 1D or 2D!')
            otherwise
                error('Wrong Input!')
        end       
        
    % 3D RVE    
    case 3
        error('3D RVE not implemented!')
    otherwise
        error('Wrong Input! Choose 2D RVE!')
end
    
    
end



