%-------------------------------------------------------------------------%
% [RVESols,MacroGPdata] = RVEPostProcess(dir_output,Mesh,step,RVEdata,
% comsolPort) is the master function that calls the appropriate RVE post
% processing function.
%
% NOTE:   Addition of new RVE models to COMSOLFE2 requires adding the
%         corresponding RVEPostProcess function to the 'switch' statement.
%
% INPUT:  dir_output     -> Output directory
%         Mesh           -> Macro-scale mesh
%         step           -> Current step number
%         RVEdata        -> MATLAB struct, must contain 'type'
%         comsolPort     -> Vector containing Port numbers associated with 
%                           the COMSOL servers
%
% OUTPUT: RVESols        -> RVE Solutions
%         MacroGPdata    -> Post-processed RVEdata (or macro-scale Gauss
%                           point data. For instance, homogenised stress,
%                           damage, etc.
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   11.01.2021
%-------------------------------------------------------------------------%

function [RVESols,MacroGPdata] = RVEPostProcess(Mesh,step,RVEdata,comsolPort)

tmp1 = cell(1,Mesh.noElements);
tmp2 = cell(1,Mesh.noElements);

switch RVEdata.type
    
    case 'ELAST1'
        
        parfor idx = 1:Mesh.noElements
            t = getCurrentTask();
            taskID=t.ID;
            try
                mphstart(comsolPort(taskID));
            catch
            end 
            import('com.comsol.model.*');
            import('com.comsol.model.util.*');
		    [tmp1{1,idx},tmp2{1,idx}] = RVEPostProcessELAST1(idx,step);
        end
        
		RVESols     = tmp1;
        MacroGPdata = tmp2;
		clear tmp1 tmp2;
    
    case 'GDDAM1'
        
        parfor idx = 1:Mesh.noElements
            t = getCurrentTask();
            taskID=t.ID;
            try
                mphstart(comsolPort(taskID));
            catch
            end 
            import('com.comsol.model.*');
            import('com.comsol.model.util.*');
		    [tmp1{1,idx},tmp2{1,idx}] = RVEPostProcessGDDAM1(idx,step);
        end
        
		RVESols     = tmp1;
        MacroGPdata = tmp2;
		clear tmp1 tmp2;
        
    case 'PFDAM1'
        
        parfor idx = 1:Mesh.noElements
            t = getCurrentTask();
            taskID=t.ID;
            try
                mphstart(comsolPort(taskID));
            catch
            end 
            import('com.comsol.model.*');
            import('com.comsol.model.util.*');
		    [tmp1{1,idx},tmp2{1,idx}] = RVEPostProcessPFDAM1(idx,step);
        end
        
		RVESols     = tmp1;
        MacroGPdata = tmp2;
		clear tmp1 tmp2;
        
    otherwise
            
        error('RVE model not implemented! Check ./src/rve/RVEInit.m')       

end

