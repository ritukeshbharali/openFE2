%-------------------------------------------------------------------------%
% [D,stress,flag] = RVESolve(idx,strain,MacroSpaceDim,action) is a master 
% function that calls the appropriate RVESolve function for macro-scale
% element 'idx', based on the RVE type selected by the user.
%
% NOTE:   Addition of new RVE models to COMSOLFE2 requires adding the
%         corresponding RVESolve function to the 'switch' statement.
%
% INPUT:  idx            -> Macro-scale element
%         strain         -> Macro-scale strain
%         MacroSpaceDim  -> Spatial Dimension of Macro-scale problem
%         action         -> 'Jac' if computation of D-tensor is required
%
% OUTPUT: D              -> Macro-scale homogenised tangent moduli
%         stress         -> Macro-scale/homogenised stress
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

function [D,stress,flag] = RVESolve(idx,strain,MacroSpaceDim,action,RVEdata)

% Selects an RVE model based on RVE 'type'
switch RVEdata.type
    
    case 'ELAST1'
        
        [D,stress,flag] = RVESolveELAST1(idx,strain,MacroSpaceDim,...
                                             action);
                                         
    case 'GDDAM1'
        
        [D,stress,flag] = RVESolveGDDAM1(idx,strain,MacroSpaceDim,...
                                             action);                                         
        
    case 'PFDAM1'
        
        [D,stress,flag] = RVESolvePFDAM1(idx,strain,MacroSpaceDim,...
                                             action);
        
    % Add RVE model name and the corresponding RVESolve function here!
    % Example:
    % case 'pfdamage2'
    %   [D,stress,flag] = RVESolvePFDAM2(idx,strain,MacroSpaceDim,...
    %                                       action);
        
    otherwise
        error('RVE model not implemented! Check ./src/rve/RVESolve.m')       
end
       
end


