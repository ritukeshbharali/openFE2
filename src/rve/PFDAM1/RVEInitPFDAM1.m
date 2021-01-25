%-------------------------------------------------------------------------%
% [D,RVEdofs] = RVEInitPFDamage1(idx,MacroSpaceDim,RVEdata) initialises a
% phase-field damage RVE with an initial slit and computes the homogenised
% tangent stiffness moduli 'D'.
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

function [D,RVEdofs] = RVEInitPFDAM1(idx,MacroSpaceDim,RVEdata)

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
        
        %    --------------------
        %    |                  |
        %    |                  |
        %    |         |        |
        %    |         |        |
        %    |                  |
        %    |                  |
        %    --------------------
        %
        %    RVE with initial crack/fracture 
        
        
        % Create a 2D RVE model
        model = ModelUtil.create('Model');
        
        % Uncomment if you wish to observe model progress in separate
        % windows. Note that one window will pop up for each process!
        % ModelUtil.showProgress(true);
        
        % Set runtime parameters (Modify as required!)
        model.param.set('Gc', '2700[N/m]', 'Critial Energy Release Rate');
        model.param.set('k', '1e-10', 'Control Numerical Sigularity');
        model.param.set('l0', '1.5e-2[mm]');
        model.param.set('hmax', 'l0/2', 'maximum element size');
        model.param.set('mu', '80.769[GPa]', 'Lame Paramter1');
        model.param.set('lanta', '131.154[GPa]', 'Lame Paramter2');
        model.param.set('eps_xx', '0');
        model.param.set('eps_yy', '0');
        model.param.set('eps_xy', '0');
        
        % Create component, geometry and mesh
        model.component.create('comp1', true);
        model.component('comp1').geom.create('geom1', 2);
        model.component('comp1').geom('geom1').lengthUnit('mm');
        model.component('comp1').geom('geom1').useConstrDim(false);
        model.component('comp1').geom('geom1').create('r1', 'Rectangle');
        model.component('comp1').geom('geom1').feature('r1').set('pos', {'-0.5[mm]' '-0.5[mm]'});
        model.component('comp1').geom('geom1').feature('r1').set('size', {'1[mm]' '1[mm]'});
        model.component('comp1').geom('geom1').create('ic1', 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature('ic1').set('table', [0.35 -0.15; 0.35 0.15 ]);
        model.component('comp1').geom('geom1').create('ic2', 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature('ic2').set('table', {'3*l0+0.35' '-0.5'; '3*l0+0.35' '0.5'});
        model.component('comp1').geom('geom1').create('ic3', 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature('ic3').set('table', {'-3*l0+0.35' '-0.5'; '-3*l0+0.35' '0.5'});
        model.component('comp1').geom('geom1').run;
        model.component('comp1').sorder('linear');
        model.component('comp1').curvedInterior(false);
        model.component('comp1').mesh.create('mesh1');
        model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
        model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
        model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
        model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([2]);
        model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
        model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
        model.component('comp1').mesh('mesh1').feature('size').set('hmax', 'hmax*8');
        model.component('comp1').mesh('mesh1').feature('size').set('hmin', 'hmax*8');
        model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('custom', 'on');
        model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmax', 'hmax');
        model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmaxactive', true);
        model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmin', 'hmax');
        model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hminactive', true);
        model.component('comp1').mesh('mesh1').run;
        
        % Some additional definitions for periodicity conditions
        model.component('comp1').selection.create('sel1', 'Explicit');
        model.component('comp1').selection('sel1').geom('geom1', 1);
        model.component('comp1').selection('sel1').set([8]);
        model.component('comp1').selection.create('sel2', 'Explicit');
        model.component('comp1').selection('sel2').geom('geom1', 1);
        model.component('comp1').selection('sel2').set([3 5 7]);
        model.component('comp1').selection.create('sel3', 'Explicit');
        model.component('comp1').selection('sel3').geom('geom1', 1);
        model.component('comp1').selection('sel3').set([2 4 6]);
        model.component('comp1').selection.create('sel4', 'Explicit');
        model.component('comp1').selection('sel4').geom('geom1', 1);
        model.component('comp1').selection('sel4').set([1]);
        model.component('comp1').selection('sel1').label('image x');
        model.component('comp1').selection('sel2').label('image y');
        model.component('comp1').selection('sel3').label('mirror y');
        model.component('comp1').selection('sel4').label('mirror x');
        model.component('comp1').cpl.create('intop1', 'Integration');
        model.component('comp1').cpl.create('linext1', 'LinearExtrusion');
        model.component('comp1').cpl.create('linext2', 'LinearExtrusion');
        model.component('comp1').cpl('intop1').selection.geom('geom1', 1);
        model.component('comp1').cpl('intop1').selection.set([3 5 7]);
        model.component('comp1').cpl('linext1').selection.geom('geom1', 1);
        model.component('comp1').cpl('linext1').selection.set([1]);
        model.component('comp1').cpl('linext2').selection.geom('geom1', 1);
        model.component('comp1').cpl('linext2').selection.set([2 4 6]);
        model.component('comp1').cpl('linext1').label('mapping x');
        model.component('comp1').cpl('linext1').set('opname', 'linextx');
        model.component('comp1').cpl('linext1').selection('srcvertex1').set([1]);
        model.component('comp1').cpl('linext1').selection('srcvertex2').set([2]);
        model.component('comp1').cpl('linext1').selection('dstvertex1').set([9]);
        model.component('comp1').cpl('linext1').selection('dstvertex2').set([10]);
        model.component('comp1').cpl('linext2').label('mapping y');
        model.component('comp1').cpl('linext2').set('opname', 'linexty');
        model.component('comp1').cpl('linext2').selection('srcvertex1').set([1]);
        model.component('comp1').cpl('linext2').selection('srcvertex2').set([9]);
        model.component('comp1').cpl('linext2').selection('dstvertex1').set([2]);
        model.component('comp1').cpl('linext2').selection('dstvertex2').set([10]);

        % Constitutive relations
        model.component('comp1').variable.create('var1');
        model.component('comp1').variable('var1').label('Constitutive Relations');
        model.component('comp1').variable('var1').set('e1_p', 'if(solid.ep1>0,solid.ep1,0)');
        model.component('comp1').variable('var1').set('e2_p', 'if(solid.ep2>0,solid.ep2,0)');
        model.component('comp1').variable('var1').set('e3_p', 'if(solid.ep3>0,solid.ep3,0)');
        model.component('comp1').variable('var1').set('tra1', 'solid.ep1+solid.ep2+solid.ep3');
        model.component('comp1').variable('var1').set('tra1_p', 'if(tra1>0,tra1,0)');
        model.component('comp1').variable('var1').set('fai_p', 'lanta*tra1_p^2/2+mu*(e1_p^2+e2_p^2+e3_p^2)');
        model.component('comp1').variable('var1').set('n1', 'if(tra1>0,(1-k)*(1-u)^2+k,1)');
        model.component('comp1').variable('var1').set('m1', 'if(solid.ep1>0,(1-k)*(1-u)^2+k,1)');
        model.component('comp1').variable('var1').set('n11', 'solid.ep1X');
        model.component('comp1').variable('var1').set('n12', 'solid.ep1Y');
        model.component('comp1').variable('var1').set('n13', 'solid.ep1Z');
        model.component('comp1').variable('var1').set('n21', 'solid.ep2X');
        model.component('comp1').variable('var1').set('n22', 'solid.ep2Y');
        model.component('comp1').variable('var1').set('n23', 'solid.ep2Z');
        model.component('comp1').variable('var1').set('n31', 'solid.ep3X');
        model.component('comp1').variable('var1').set('n32', 'solid.ep3Y');
        model.component('comp1').variable('var1').set('n33', 'solid.ep3Z');
        model.component('comp1').variable('var1').set('d11_p', 'if(solid.ep1>0,1,0)');
        model.component('comp1').variable('var1').set('d22_p', 'if(solid.ep2>0,1,0)');
        model.component('comp1').variable('var1').set('d33_p', 'if(solid.ep3>0,1,0)');
        model.component('comp1').variable('var1').set('d11_n', 'if(solid.ep1<=0,1,0)');
        model.component('comp1').variable('var1').set('d22_n', 'if(solid.ep2<=0,1,0)');
        model.component('comp1').variable('var1').set('d33_n', 'if(solid.ep3<=0,1,0)');
        model.component('comp1').variable('var1').set('l1', 'if(solid.ep1>solid.ep2,solid.ep1,solid.ep1+1e-9)');
        model.component('comp1').variable('var1').set('l2', 'solid.ep2');
        model.component('comp1').variable('var1').set('l3', 'if(solid.ep3<solid.ep2,solid.ep3,solid.ep3-1e-9)');
        model.component('comp1').variable('var1').set('g1_p', 'if(solid.ep1>solid.ep2,d11_p*solid.ep1,d11_p*solid.ep1+1e-9)');
        model.component('comp1').variable('var1').set('g2_p', 'd22_p*solid.ep2');
        model.component('comp1').variable('var1').set('g3_p', 'if(solid.ep3<solid.ep2,d33_p*solid.ep3,d33_p*solid.ep3-1e-9)');
        model.component('comp1').variable('var1').set('g1_n', 'if(solid.ep1>solid.ep2,d11_n*solid.ep1,d11_n*solid.ep1+1e-9)');
        model.component('comp1').variable('var1').set('g2_n', 'd22_n*solid.ep2');
        model.component('comp1').variable('var1').set('g3_n', 'if(solid.ep3<solid.ep2,d33_n*solid.ep3,d33_n*solid.ep3-1e-9)');
        model.component('comp1').variable('var1').set('g1111_p', 'd11_p*n11*n11*n11*n11+d22_p*n21*n21*n21*n21+d33_p*n31*n31*n31*n31+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n11*n21+n21*n11)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n11*n31+n31*n11)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n21*n11+n11*n21)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n21*n31+n31*n21)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n31*n11+n11*n31)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n31*n21+n21*n31)');
        model.component('comp1').variable('var1').set('g1122_p', 'd11_p*n11*n11*n12*n12+d22_p*n21*n21*n22*n22+d33_p*n31*n31*n32*n32+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n12*n22+n22*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n12*n32+n32*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n22*n12+n12*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n22*n32+n32*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n32*n12+n12*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n32*n22+n22*n32)');
        model.component('comp1').variable('var1').set('g1133_p', 'd11_p*n11*n11*n13*n13+d22_p*n21*n21*n23*n23+d33_p*n31*n31*n33*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n13*n23+n23*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n13*n33+n33*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n23*n13+n13*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n23*n33+n33*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n33*n13+n13*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n33*n23+n23*n33)');
        model.component('comp1').variable('var1').set('g1112_p', 'd11_p*n11*n11*n11*n12+d22_p*n21*n21*n21*n22+d33_p*n31*n31*n31*n32+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n11*n22+n21*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n11*n32+n31*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n21*n12+n11*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n21*n32+n31*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n31*n12+n11*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n31*n22+n21*n32)');
        model.component('comp1').variable('var1').set('g1123_p', 'd11_p*n11*n11*n12*n13+d22_p*n21*n21*n22*n23+d33_p*n31*n31*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g1113_p', 'd11_p*n11*n11*n11*n13+d22_p*n21*n21*n21*n23+d33_p*n31*n31*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n21*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n31*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n11*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n31*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n11*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n21*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g2222_p', 'd11_p*n12*n12*n12*n12+d22_p*n22*n22*n22*n22+d33_p*n32*n32*n32*n32+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n12*n22+n22*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n12*n32+n32*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n22*n12+n12*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n22*n32+n32*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n32*n12+n12*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n32*n22+n22*n32)');
        model.component('comp1').variable('var1').set('g2233_p', 'd11_p*n12*n12*n13*n13+d22_p*n22*n22*n23*n23+d33_p*n32*n32*n33*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n13*n23+n23*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n13*n33+n33*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n23*n13+n13*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n23*n33+n33*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n33*n13+n13*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n33*n23+n23*n33)');
        model.component('comp1').variable('var1').set('g2212_p', 'd11_p*n12*n12*n11*n12+d22_p*n22*n22*n21*n22+d33_p*n32*n32*n31*n32+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n11*n22+n21*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n11*n32+n31*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n21*n12+n11*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n21*n32+n31*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n31*n12+n11*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n31*n22+n21*n32)');
        model.component('comp1').variable('var1').set('g2223_p', 'd11_p*n12*n12*n12*n13+d22_p*n22*n22*n22*n23+d33_p*n32*n32*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g2213_p', 'd11_p*n12*n12*n11*n13+d22_p*n22*n22*n21*n23+d33_p*n32*n32*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n22*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n32*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n12*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n32*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n12*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n22*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g3333_p', 'd11_p*n13*n13*n13*n13+d22_p*n23*n23*n23*n23+d33_p*n33*n33*n33*n33+0.5*(g1_p-g2_p)/(l1-l2)*n13*n23*(n13*n23+n23*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n13*n33*(n13*n33+n33*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n23*n13*(n23*n13+n13*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n23*n33*(n23*n33+n33*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n33*n13*(n33*n13+n13*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n33*n23*(n33*n23+n23*n33)');
        model.component('comp1').variable('var1').set('g3312_p', 'd11_p*n13*n13*n11*n12+d22_p*n23*n23*n21*n22+d33_p*n33*n33*n31*n32+0.5*(g1_p-g2_p)/(l1-l2)*n13*n23*(n11*n22+n21*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n13*n33*(n11*n32+n31*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n23*n13*(n21*n12+n11*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n23*n33*(n21*n32+n31*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n33*n13*(n31*n12+n11*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n33*n23*(n31*n22+n21*n32)');
        model.component('comp1').variable('var1').set('g3323_p', 'd11_p*n13*n13*n12*n13+d22_p*n23*n23*n22*n23+d33_p*n33*n33*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n13*n23*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n13*n33*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n23*n13*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n23*n33*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n33*n13*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n33*n23*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g3313_p', 'd11_p*n13*n13*n11*n13+d22_p*n23*n23*n21*n23+d33_p*n33*n33*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n13*n23*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n13*n33*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n23*n13*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n23*n33*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n33*n13*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n33*n23*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g1212_p', 'd11_p*n11*n12*n11*n12+d22_p*n21*n22*n21*n22+d33_p*n31*n32*n31*n32+0.5*(g1_p-g2_p)/(l1-l2)*n11*n22*(n11*n22+n21*n12)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n32*(n11*n32+n31*n12)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n12*(n21*n12+n11*n22)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n32*(n21*n32+n31*n22)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n12*(n31*n12+n11*n32)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n22*(n31*n22+n21*n32)');
        model.component('comp1').variable('var1').set('g1223_p', 'd11_p*n11*n12*n12*n13+d22_p*n21*n22*n22*n23+d33_p*n31*n32*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n22*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n32*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n12*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n32*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n12*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n22*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g1213_p', 'd11_p*n11*n12*n11*n13+d22_p*n21*n22*n21*n23+d33_p*n31*n32*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n22*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n32*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n12*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n32*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n12*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n22*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g2323_p', 'd11_p*n12*n13*n12*n13+d22_p*n22*n23*n22*n23+d33_p*n32*n33*n32*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n23*(n12*n23+n22*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n33*(n12*n33+n32*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n13*(n22*n13+n12*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n33*(n22*n33+n32*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n13*(n32*n13+n12*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n23*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g2313_p', 'd11_p*n12*n13*n11*n13+d22_p*n22*n23*n21*n23+d33_p*n32*n33*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n12*n23*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n12*n33*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n22*n13*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n22*n33*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n32*n13*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n32*n23*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g1313_p', 'd11_p*n11*n13*n11*n13+d22_p*n21*n23*n21*n23+d33_p*n31*n33*n31*n33+0.5*(g1_p-g2_p)/(l1-l2)*n11*n23*(n11*n23+n21*n13)+0.5*(g1_p-g3_p)/(l1-l3)*n11*n33*(n11*n33+n31*n13)+0.5*(g2_p-g1_p)/(l2-l1)*n21*n13*(n21*n13+n11*n23)+0.5*(g2_p-g3_p)/(l2-l3)*n21*n33*(n21*n33+n31*n23)+0.5*(g3_p-g1_p)/(l3-l1)*n31*n13*(n31*n13+n11*n33)+0.5*(g3_p-g2_p)/(l3-l2)*n31*n23*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g1111_n', 'd11_n*n11*n11*n11*n11+d22_n*n21*n21*n21*n21+d33_n*n31*n31*n31*n31+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n11*n21+n21*n11)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n11*n31+n31*n11)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n21*n11+n11*n21)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n21*n31+n31*n21)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n31*n11+n11*n31)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n31*n21+n21*n31)');
        model.component('comp1').variable('var1').set('g1122_n', 'd11_n*n11*n11*n12*n12+d22_n*n21*n21*n22*n22+d33_n*n31*n31*n32*n32+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n12*n22+n22*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n12*n32+n32*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n22*n12+n12*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n22*n32+n32*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n32*n12+n12*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n32*n22+n22*n32)');
        model.component('comp1').variable('var1').set('g1133_n', 'd11_n*n11*n11*n13*n13+d22_n*n21*n21*n23*n23+d33_n*n31*n31*n33*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n13*n23+n23*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n13*n33+n33*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n23*n13+n13*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n23*n33+n33*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n33*n13+n13*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n33*n23+n23*n33)');
        model.component('comp1').variable('var1').set('g1112_n', 'd11_n*n11*n11*n11*n12+d22_n*n21*n21*n21*n22+d33_n*n31*n31*n31*n32+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n11*n22+n21*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n11*n32+n31*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n21*n12+n11*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n21*n32+n31*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n31*n12+n11*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n31*n22+n21*n32)');
        model.component('comp1').variable('var1').set('g1123_n', 'd11_n*n11*n11*n12*n13+d22_n*n21*n21*n22*n23+d33_n*n31*n31*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g1113_n', 'd11_n*n11*n11*n11*n13+d22_n*n21*n21*n21*n23+d33_n*n31*n31*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n21*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n31*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n11*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n31*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n11*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n21*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g2222_n', 'd11_n*n12*n12*n12*n12+d22_n*n22*n22*n22*n22+d33_n*n32*n32*n32*n32+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n12*n22+n22*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n12*n32+n32*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n22*n12+n12*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n22*n32+n32*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n32*n12+n12*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n32*n22+n22*n32)');
        model.component('comp1').variable('var1').set('g2233_n', 'd11_n*n12*n12*n13*n13+d22_n*n22*n22*n23*n23+d33_n*n32*n32*n33*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n13*n23+n23*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n13*n33+n33*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n23*n13+n13*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n23*n33+n33*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n33*n13+n13*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n33*n23+n23*n33)');
        model.component('comp1').variable('var1').set('g2212_n', 'd11_n*n12*n12*n11*n12+d22_n*n22*n22*n21*n22+d33_n*n32*n32*n31*n32+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n11*n22+n21*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n11*n32+n31*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n21*n12+n11*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n21*n32+n31*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n31*n12+n11*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n31*n22+n21*n32)');
        model.component('comp1').variable('var1').set('g2223_n', 'd11_n*n12*n12*n12*n13+d22_n*n22*n22*n22*n23+d33_n*n32*n32*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g2213_n', 'd11_n*n12*n12*n11*n13+d22_n*n22*n22*n21*n23+d33_n*n32*n32*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n22*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n32*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n12*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n32*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n12*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n22*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g3333_n', 'd11_n*n13*n13*n13*n13+d22_n*n23*n23*n23*n23+d33_n*n33*n33*n33*n33+0.5*(g1_n-g2_n)/(l1-l2)*n13*n23*(n13*n23+n23*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n13*n33*(n13*n33+n33*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n23*n13*(n23*n13+n13*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n23*n33*(n23*n33+n33*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n33*n13*(n33*n13+n13*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n33*n23*(n33*n23+n23*n33)');
        model.component('comp1').variable('var1').set('g3312_n', 'd11_n*n13*n13*n11*n12+d22_n*n23*n23*n21*n22+d33_n*n33*n33*n31*n32+0.5*(g1_n-g2_n)/(l1-l2)*n13*n23*(n11*n22+n21*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n13*n33*(n11*n32+n31*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n23*n13*(n21*n12+n11*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n23*n33*(n21*n32+n31*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n33*n13*(n31*n12+n11*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n33*n23*(n31*n22+n21*n32)');
        model.component('comp1').variable('var1').set('g3323_n', 'd11_n*n13*n13*n12*n13+d22_n*n23*n23*n22*n23+d33_n*n33*n33*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n13*n23*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n13*n33*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n23*n13*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n23*n33*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n33*n13*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n33*n23*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g3313_n', 'd11_n*n13*n13*n11*n13+d22_n*n23*n23*n21*n23+d33_n*n33*n33*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n13*n23*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n13*n33*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n23*n13*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n23*n33*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n33*n13*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n33*n23*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g1212_n', 'd11_n*n11*n12*n11*n12+d22_n*n21*n22*n21*n22+d33_n*n31*n32*n31*n32+0.5*(g1_n-g2_n)/(l1-l2)*n11*n22*(n11*n22+n21*n12)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n32*(n11*n32+n31*n12)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n12*(n21*n12+n11*n22)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n32*(n21*n32+n31*n22)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n12*(n31*n12+n11*n32)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n22*(n31*n22+n21*n32)');
        model.component('comp1').variable('var1').set('g1223_n', 'd11_n*n11*n12*n12*n13+d22_n*n21*n22*n22*n23+d33_n*n31*n32*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n22*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n32*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n12*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n32*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n12*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n22*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g1213_n', 'd11_n*n11*n12*n11*n13+d22_n*n21*n22*n21*n23+d33_n*n31*n32*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n22*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n32*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n12*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n32*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n12*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n22*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g2323_n', 'd11_n*n12*n13*n12*n13+d22_n*n22*n23*n22*n23+d33_n*n32*n33*n32*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n23*(n12*n23+n22*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n33*(n12*n33+n32*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n13*(n22*n13+n12*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n33*(n22*n33+n32*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n13*(n32*n13+n12*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n23*(n32*n23+n22*n33)');
        model.component('comp1').variable('var1').set('g2313_n', 'd11_n*n12*n13*n11*n13+d22_n*n22*n23*n21*n23+d33_n*n32*n33*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n12*n23*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n12*n33*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n22*n13*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n22*n33*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n32*n13*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n32*n23*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('g1313_n', 'd11_n*n11*n13*n11*n13+d22_n*n21*n23*n21*n23+d33_n*n31*n33*n31*n33+0.5*(g1_n-g2_n)/(l1-l2)*n11*n23*(n11*n23+n21*n13)+0.5*(g1_n-g3_n)/(l1-l3)*n11*n33*(n11*n33+n31*n13)+0.5*(g2_n-g1_n)/(l2-l1)*n21*n13*(n21*n13+n11*n23)+0.5*(g2_n-g3_n)/(l2-l3)*n21*n33*(n21*n33+n31*n23)+0.5*(g3_n-g1_n)/(l3-l1)*n31*n13*(n31*n13+n11*n33)+0.5*(g3_n-g2_n)/(l3-l2)*n31*n23*(n31*n23+n21*n33)');
        model.component('comp1').variable('var1').set('cc', '(1-k)*(1-u)^2+k');
        model.component('comp1').variable('var1').set('dd', 'if(tra1>0,(1-k)*(1-u)^2+k,1)');
        model.component('comp1').variable('var1').set('g1111', '2*mu*(cc*g1111_p+g1111_n)+lanta*dd');
        model.component('comp1').variable('var1').set('g1122', '2*mu*(cc*g1122_p+g1122_n)+lanta*dd');
        model.component('comp1').variable('var1').set('g1133', '2*mu*(cc*g1133_p+g1133_n)+lanta*dd');
        model.component('comp1').variable('var1').set('g1112', '2*mu*(cc*g1112_p+g1112_n)');
        model.component('comp1').variable('var1').set('g1123', '2*mu*(cc*g1123_p+g1123_n)');
        model.component('comp1').variable('var1').set('g1113', '2*mu*(cc*g1113_p+g1113_n)');
        model.component('comp1').variable('var1').set('g2222', '2*mu*(cc*g2222_p+g2222_n)+lanta*dd');
        model.component('comp1').variable('var1').set('g2233', '2*mu*(cc*g2233_p+g2233_n)+lanta*dd');
        model.component('comp1').variable('var1').set('g2212', '2*mu*(cc*g2212_p+g2212_n)');
        model.component('comp1').variable('var1').set('g2223', '2*mu*(cc*g2223_p+g2223_n)');
        model.component('comp1').variable('var1').set('g2213', '2*mu*(cc*g2213_p+g2213_n)');
        model.component('comp1').variable('var1').set('g3333', '2*mu*(cc*g3333_p+g3333_n)+lanta*dd');
        model.component('comp1').variable('var1').set('g3312', '2*mu*(cc*g3312_p+g3312_n)');
        model.component('comp1').variable('var1').set('g3323', '2*mu*(cc*g3323_p+g3323_n)');
        model.component('comp1').variable('var1').set('g3313', '2*mu*(cc*g3313_p+g3313_n)');
        model.component('comp1').variable('var1').set('g1212', '2*mu*(cc*g1212_p+g1212_n)');
        model.component('comp1').variable('var1').set('g1223', '2*mu*(cc*g1223_p+g1223_n)');
        model.component('comp1').variable('var1').set('g1213', '2*mu*(cc*g1213_p+g1213_n)');
        model.component('comp1').variable('var1').set('g2323', '2*mu*(cc*g2323_p+g2323_n)');
        model.component('comp1').variable('var1').set('g2313', '2*mu*(cc*g2313_p+g2313_n)');
        model.component('comp1').variable('var1').set('g1313', '2*mu*(cc*g1313_p+g1313_n)');
        model.component('comp1').variable('var1').set('rex', 'intop1(solid.sxy)');
        model.component('comp1').variable('var1').selection.geom('geom1', 2);
        model.component('comp1').variable('var1').selection.all;

        % Periodicity in x-direction
        model.component('comp1').variable.create('var2');
        model.component('comp1').variable('var2').label('on positive x');
        model.component('comp1').variable('var2').set('dx', 'x-linextx(x)');
        model.component('comp1').variable('var2').set('dy', 'y-linextx(y)');
        model.component('comp1').variable('var2').set('du2', 'u2-linextx(u2)');
        model.component('comp1').variable('var2').set('dv2', 'v2-linextx(v2)');
        model.component('comp1').variable('var2').set('du', 'u-linextx(u)');
        model.component('comp1').variable('var2').selection.geom('geom1', 1);
        model.component('comp1').variable('var2').selection.set([8]);

        % Periodicity in y-direction
        model.component('comp1').variable.create('var3');
        model.component('comp1').variable('var3').label('on positive y');
        model.component('comp1').variable('var3').set('dx', 'x-linexty(x)');
        model.component('comp1').variable('var3').set('dy', 'y-linexty(y)');
        model.component('comp1').variable('var3').set('du2', 'u2-linexty(u2)');
        model.component('comp1').variable('var3').set('dv2', 'v2-linexty(v2)');
        model.component('comp1').variable('var3').set('du', 'u-linexty(u)');
        model.component('comp1').variable('var3').selection.geom('geom1', 1);
        model.component('comp1').variable('var3').selection.set([3 5 7]);

        % Solid Mechanics module (Momentum balance equation)
        model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');
        model.component('comp1').physics('solid').field('displacement').field('u2');
        model.component('comp1').physics('solid').field('displacement').component({'u2' 'v2' 'w2'});
        model.component('comp1').physics('solid').create('constr1', 'PointwiseConstraint', 1);
        model.component('comp1').physics('solid').feature('constr1').selection.set([3 5 7 8]);
        model.component('comp1').physics('solid').create('constr2', 'PointwiseConstraint', 1);
        model.component('comp1').physics('solid').feature('constr2').selection.set([3 5 7 8]);
        model.component('comp1').physics('solid').create('fix1', 'Fixed', 0);
        model.component('comp1').physics('solid').feature('fix1').selection.set([1]);
        model.component('comp1').physics('solid').prop('ShapeProperty').set('order_displacement', 1);
        model.component('comp1').physics('solid').prop('StructuralTransientBehavior').set('StructuralTransientBehavior', 'Quasistatic');
        model.component('comp1').physics('solid').feature('lemm1').set('SolidModel', 'Anisotropic');
        model.component('comp1').physics('solid').feature('lemm1').set('D_mat', 'userdef');
        model.component('comp1').physics('solid').feature('lemm1').set('D', {'g1111';  ...
        'g1122';  ...
        'g1133';  ...
        'g1112';  ...
        'g1123';  ...
        'g1113';  ...
        'g1122';  ...
        'g2222';  ...
        'g2233';  ...
        'g2212';  ...
        'g2223';  ...
        'g2213';  ...
        'g1133';  ...
        'g2233';  ...
        'g3333';  ...
        'g3312';  ...
        'g3323';  ...
        'g3313';  ...
        'g1112';  ...
        'g2212';  ...
        'g3312';  ...
        'g1212';  ...
        'g1223';  ...
        'g1213';  ...
        'g1123';  ...
        'g2223';  ...
        'g3323';  ...
        'g1223';  ...
        'g2323';  ...
        'g2313';  ...
        'g1113';  ...
        'g2213';  ...
        'g3313';  ...
        'g1213';  ...
        'g2313';  ...
        'g1313'});
        model.component('comp1').physics('solid').feature('lemm1').set('rho_mat', 'userdef');
        model.component('comp1').physics('solid').feature('lemm1').set('rho', 0);
        model.component('comp1').physics('solid').feature('lemm1').set('minput_temperature_src', 'userdef');
        model.component('comp1').physics('solid').feature('constr1').set('constraintExpression', '-du2+eps_xx*dx+eps_xy*dy');
        model.component('comp1').physics('solid').feature('constr1').set('shapeFunctionType', 'shgp');
        model.component('comp1').physics('solid').feature('constr1').set('constraintMethod', 'elemental');
        model.component('comp1').physics('solid').feature('constr1').set('order', 4);
        model.component('comp1').physics('solid').feature('constr1').set('frame', 'spatial');
        model.component('comp1').physics('solid').feature('constr2').set('constraintExpression', '-dv2+eps_xy*dx+eps_yy*dy');
        model.component('comp1').physics('solid').feature('constr2').set('shapeFunctionType', 'shgp');
        model.component('comp1').physics('solid').feature('constr2').set('order', 4);
        model.component('comp1').physics('solid').feature('constr2').set('constraintMethod', 'elemental');
        model.component('comp1').physics('solid').feature('constr2').set('frame', 'spatial');

        % Weak Form module (Phase-field evolution equation)
        model.component('comp1').physics.create('w', 'WeakFormPDE', 'geom1');
        model.component('comp1').physics('w').create('constr1', 'PointwiseConstraint', 1);
        model.component('comp1').physics('w').feature('constr1').selection.set([3 5 7 8]);
        model.component('comp1').physics('w').create('dir1', 'DirichletBoundary', 1);
        model.component('comp1').physics('w').feature('dir1').selection.set([10]);
        model.component('comp1').physics('w').prop('ShapeProperty').set('order', 1);
        model.component('comp1').physics('w').prop('Units').set('CustomSourceTermUnit', 'N/mm^2');
        model.component('comp1').physics('w').feature('wfeq1').set('weak', '-Gc*l0*(test(ux)*ux+test(uy)*uy) + test(u)*2*(1-u)*H - test(u)*Gc*u/l0');
        model.component('comp1').physics('w').feature('constr1').set('constraintExpression', '-du');
        model.component('comp1').physics('w').feature('constr1').set('constraintMethod', 'elemental');
        model.component('comp1').physics('w').feature('constr1').set('shapeFunctionType', 'shgp');
        model.component('comp1').physics('w').feature('constr1').set('order', 4);
        model.component('comp1').physics('w').feature('dir1').set('r', 1);
        model.component('comp1').physics('w').feature('dir1').set('constraintType', 'symmetricConstraint');
        model.component('comp1').physics('w').feature('dir1').set('constraintMethod', 'nodal');

        % Domain ODE module (Update History-variable [Miehe,2010])
        model.component('comp1').physics.create('dode', 'DomainODE', 'geom1');
        model.component('comp1').physics('dode').field('dimensionless').field('H');
        model.component('comp1').physics('dode').field('dimensionless').component({'H'});
        model.component('comp1').physics('dode').prop('Units').set('DependentVariableQuantity', 'none');
        model.component('comp1').physics('dode').prop('Units').set('CustomDependentVariableUnit', 'N/mm^2');
        model.component('comp1').physics('dode').prop('ShapeProperty').set('shapeFunctionType', 'shgp');
        model.component('comp1').physics('dode').prop('ShapeProperty').set('order', 4);
        model.component('comp1').physics('dode').prop('Units').set('CustomSourceTermUnit', 'N/mm^2');
        model.component('comp1').physics('dode').feature('dode1').set('f', 'H-nojac(if(fai_p>H,fai_p,H))');
        model.component('comp1').physics('dode').feature('dode1').set('da', 0);

        % Define solver
        model.study.create('std2');
        model.study('std2').create('time', 'Transient');
        model.sol.create('sol1');
        model.sol('sol1').study('std2');
        model.sol('sol1').attach('std2');
        model.sol('sol1').create('st1', 'StudyStep');
        model.sol('sol1').create('v1', 'Variables');
        model.sol('sol1').create('t1', 'Time');
        model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
        model.sol('sol1').feature('t1').create('se1', 'Segregated');
        model.sol('sol1').feature('t1').create('ps1', 'PreviousSolution');
        model.sol('sol1').feature('t1').create('st1', 'StopCondition');
        model.sol('sol1').feature('t1').feature('se1').create('ss1', 'SegregatedStep');
        model.sol('sol1').feature('t1').feature('se1').create('ss2', 'SegregatedStep');
        model.sol('sol1').feature('t1').feature('se1').create('ll1', 'LowerLimit');
        model.sol('sol1').feature('t1').feature('se1').create('ul1', 'UpperLimit');
        model.sol('sol1').feature('t1').feature.remove('fcDef');
        model.study('std2').feature('time').set('tlist', 'range(0,0.5,1)');
        model.study('std2').feature('time').set('usertol', true);
        model.study('std2').feature('time').set('rtol', '1e-3');
        model.study('std2').feature('time').set('useinitsol', true);
        model.study('std2').feature('time').set('initmethod', 'sol');    
        model.sol('sol1').attach('std2');
        model.sol('sol1').feature('v1').set('initmethod', 'sol');
        model.sol('sol1').feature('v1').set('clist', {'range(0,0.5,1)' '0.001[s]'});
        model.sol('sol1').feature('v1').feature('comp1_u2').set('scalemethod', 'manual');
        model.sol('sol1').feature('v1').feature('comp1_u2').set('scaleval', '1e-2*0.0014142135623730952');
        model.sol('sol1').feature('t1').set('control', 'user');
        model.sol('sol1').feature('t1').set('tlist', 'range(0,0.5,1)');
        model.sol('sol1').feature('t1').set('rtol', '1e-3');
        model.sol('sol1').feature('t1').set('timemethod', 'genalpha');
        model.sol('sol1').feature('t1').set('rhoinf', 0.85);
        model.sol('sol1').feature('t1').set('estrat', 'exclude');
        model.sol('sol1').feature('t1').set('maxstepconstraintgenalpha', 'const');
        model.sol('sol1').feature('t1').set('maxstepgenalpha', 1);
        model.sol('sol1').feature('t1').set('rescaleafterinitbw', true);
        model.sol('sol1').feature('t1').feature('dDef').set('ooc', false);
        model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
        model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 100);
        model.sol('sol1').feature('t1').feature('fc1').set('stabacc', 'aacc');
        model.sol('sol1').feature('t1').feature('fc1').set('aaccdim', 50);
        model.sol('sol1').feature('t1').feature('se1').set('maxsegiter', 500);
        model.sol('sol1').feature('t1').feature('se1').set('segstabacc', 'segaacc');
        model.sol('sol1').feature('t1').feature('se1').set('segaaccdim', 100);
        model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_u2'});
        model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('segvar', {'comp1_H'});
        model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('segvar', {'comp1_u'});
        model.sol('sol1').feature('t1').feature('ps1').set('prevcomp', {'comp1_H'});
        model.sol('sol1').feature('t1').feature('st1').set('stopcondterminateon', {'true'});
        model.sol('sol1').feature('t1').feature('st1').set('stopcondActive', {'on'});
        model.sol('sol1').feature('t1').feature('st1').set('stopconddesc', {'Stop expression 1'});
        model.sol('sol1').feature('t1').feature('st1').set('stopcondarr', {'timestep<1e-6[s]'});
        
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
        
        % Create post-processing step to compute the failure-zone averaged
        % phase-field damage (TO-DO)
        
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



