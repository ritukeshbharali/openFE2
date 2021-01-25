%-------------------------------------------------------------------------%
% [RVESols,MacroGPdata] = RVEPostProcessGDDAM1(idx,step,dir_output)
% performs a post-processing step. The function returns the RVE solutions
% 'RVESols' for all RVEs along with user-defined quantities in
% 'MacroGPdata'. Examples of the later maybe homogenised stress, damage
% etc. Additionally, the user may also print RVE image files w.r.t. a
% quantity of interest, say phase-field. 
%
% INPUT:  idx            -> Macro-scale Element
%         step           -> Current step number
%         dir_output     -> Output directory
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%
% TO-DO: Add another input parameter to control post processing requests
%         
% DATE:   21.01.2021
%-------------------------------------------------------------------------% 

function [RVESols,MacroGPdata] = RVEPostProcessGDDAM1(idx,step)

% Load the RVE corresponding to the element 'idx'
RVEname = sprintf('RVE%d.mph', idx);
folder = './output/RVEs';
model = mphload(fullfile(folder,RVEname));

dir_postproc = './output/images';
if ~exist(dir_postproc, 'dir')
    mkdir(dir_postproc)
end

% Image names
RVEfigsmises  = sprintf('Smises_Step%d_RVE%d', step , idx);
RVEfigsxx     = sprintf('Sxx_Step%d_RVE%d', step , idx);
RVEfigsyy     = sprintf('Syy_Step%d_RVE%d', step , idx);
RVEfigsxy     = sprintf('Sxy_Step%d_RVE%d', step , idx);

try   
    model.result('pg1').feature('surf1').set('expr','solid.sx');
    model.result('pg1').setIndex('looplevel', 2, 0);
    model.result('pg1').run;
    model.result('pg1').feature('surf1').set('resolution', 'normal');
    model.result('pg1').feature('surf1').set('colortable', 'Rainbow');
    model.result('pg1').feature('surf1').set('coloring', 'colortable');
    model.result('pg1').feature('surf1').set('rangecoloractive', false);
    f1 = figure('visible','off');
    mphplot(model,'pg1','rangenum',2);
    saveas(f1,fullfile(dir_postproc,RVEfigsxx),'png');
    
    model.result('pg1').feature('surf1').set('expr','solid.sy');
    model.result('pg1').setIndex('looplevel', 2, 0);
    model.result('pg1').run;
    model.result('pg1').feature('surf1').set('resolution', 'normal');
    model.result('pg1').feature('surf1').set('colortable', 'Rainbow');
    model.result('pg1').feature('surf1').set('coloring', 'colortable');
    model.result('pg1').feature('surf1').set('rangecoloractive', false);
    f2 = figure('visible','off');
    mphplot(model,'pg1','rangenum',2);
    saveas(f2,fullfile(dir_postproc,RVEfigsyy),'png');
    
    model.result('pg1').feature('surf1').set('expr','solid.sxy');
    model.result('pg1').setIndex('looplevel', 2, 0);
    model.result('pg1').run;
    model.result('pg1').feature('surf1').set('resolution', 'normal');
    model.result('pg1').feature('surf1').set('colortable', 'Rainbow');
    model.result('pg1').feature('surf1').set('coloring', 'colortable');
    model.result('pg1').feature('surf1').set('rangecoloractive', false);
    f3 = figure('visible','off');
    mphplot(model,'pg1','rangenum',2);
    saveas(f3,fullfile(dir_postproc,RVEfigsxy),'png');
    
    % Phase-field
    model.result('pg1').feature('surf1').set('expr','solid.dmg');
    model.result('pg1').setIndex('looplevel', 2, 0);
    model.result('pg1').run;
    model.result('pg1').feature('surf1').set('resolution', 'normal');
    model.result('pg1').feature('surf1').set('colortable', 'Rainbow');
    model.result('pg1').feature('surf1').set('coloring', 'colortable');
    model.result('pg1').feature('surf1').set('rangecoloractive', false);
    f4 = figure('visible','off');
    mphplot(model,'pg1','rangenum',2);
    saveas(f4,fullfile(dir_postproc,RVEfigsmises),'png');
    
catch
end

% Extract the 'last' RVE solution
solnums = model.sol('sol1').getPVals();
RVESols = model.sol('sol1').getU(length(solnums));

% Compute additional data (say failure zone averaged damage)
MacroGPdata = [];

% Save the RVE
mphsave(model,fullfile(folder,RVEname));    

end

