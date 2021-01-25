%-------------------------------------------------------------------------%
% [D,RVEdofs] = RVEInit(idx,MacroSpaceDim,RVEdata) is a master function
% that calls the appropriate RVE initialisation function for macro-scale
% element 'idx', based on the RVE type selected by the user.
%
% NOTE:   Addition of new RVE models to COMSOLFE2 requires adding the
%         corresponding RVEInit function to the 'switch' statement.
%
% INPUT:  idx            -> Macro-scale Element
%         MacroSpaceDim  -> Spatial Dimension of Macro-scale problem
%         RVEdata        -> MATLAB struct must contain RVE 'type'
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
%-------------------------------------------------------------------------%

function [D,RVEdofs] = RVEInit(idx,MacroSpaceDim,RVEdata)

% Selects an RVE model based on RVE 'type'
switch RVEdata.type
    
    case 'ELAST1' 
        
        [D,RVEdofs] = RVEInitELAST1(idx,MacroSpaceDim,RVEdata);
        
    case 'GDDAM1' 
        
        [D,RVEdofs] = RVEInitGDDAM1(idx,MacroSpaceDim,RVEdata);
        
    case 'PFDAM1'
        
        [D,RVEdofs] = RVEInitPFDAM1(idx,MacroSpaceDim,RVEdata);
        
    % Add RVE model name and the corresponding RVEInit function here!
    % Example:
    % case 'pfdamage2'
    %   [D,RVEdofs] = RVEInitPFDAM2(idx,MacroSpaceDim,RVEdata);
        
    otherwise
        error('RVE model not implemented! Check ./src/rve/RVEInit.m')       
end
       
end


