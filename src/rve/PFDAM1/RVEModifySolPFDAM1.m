%-------------------------------------------------------------------------%
% RVEModifySolPFDamage1(idx,sol,action) performs modifications to
% the RVE solution vector for the macro-scale element 'idx'. Apart from
% 'idx', the function requires an RVE solution vector 'sol' and 'action'
% that indicates which type of modification is required.
%
% INPUT:  idx            -> Macro-scale Element
%         sol            -> RVE solution vector
%         action         -> 'zero' sets RVE solution to zero
%                        -> 'history' sets the phase-field history variable
%                           from 'sol'
%                        -> 'revert' sets the entire RVE solution vector to
%                           'sol'. Useful, in case RVE problems fails and
%                           the macro-problem and the RVES requires 
%                           reverting to an earlier converged state.
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   06.01.2021
%-------------------------------------------------------------------------%

function RVEModifySolPFDAM1(elem_idx,sol,action)

% Load the RVE corresponding to the element 'elem_idx'
RVEname = sprintf('RVE%d.mph', elem_idx);
folder = './output/RVEs';
model = mphload(fullfile(folder,RVEname));

switch action
    case 'zero'
        tmp = mphgetu(model,'soltag','sol1');
        tmp = zeros(length(tmp),1);
        PVal = model.sol('sol1').getPVals();
        for i = 1:length(PVal)
        model.sol('sol1').setU(tmp);
        model.sol('sol1').setPVals( PVal(i) );       % Need to check this!
        model.sol('sol1').createSolution();
        end
        
    case 'history'
        % Get current solution
        tmp = mphgetu(model,'soltag','sol1');
        % Get more mesh info
        info = mphxmeshinfo(model);
        % Get index of H
        idxH = find(strcmp(info.dofs.dofnames,'comp1.H'))-1;
        % Get H dofs
        Hdofs = find(info.dofs.nameinds==idxH);
        % Set History variable dofs to old timestep values
        tmp(Hdofs) = sol(Hdofs);
        % Push modified solution ('sol1') to COMSOL
        model.sol('sol1').setU(tmp);
        model.sol('sol1').setPVals( 2 );       % Need to check this!
        model.sol('sol1').createSolution();
        
    case 'revert'
        model.sol('sol1').setU(sol);
        model.sol('sol1').setPVals( 2 );    % Requires TESTING!
        model.sol('sol1').createSolution();
    otherwise
        error('Wrong option!')

end

mphsave(model,fullfile(folder,RVEname));

end

