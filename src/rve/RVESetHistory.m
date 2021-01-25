%-------------------------------------------------------------------------%
% RVERevert(Mesh,RVEtype,micro_tmp,comsolPort) is a master function that
% calls the appropriate function to modify the history variable to a
% previously converged state. This is useful in nonlinear failure analysis,
% for instance, phase-field damage modelling.
%
% NOTE:   Addition of new RVE models to COMSOLFE2 requires adding the
%         corresponding RVEModifySol function to the 'switch' statement.
%
% INPUT:  Mesh           -> Macro-scale mesh
%         RVEtype        -> RVE type
%         micro_tmp      -> old converged RVE solution
%         comsolPort     -> Vector containing Port numbers associated with 
%                           the COMSOL servers
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   11.01.2021
%-------------------------------------------------------------------------%

function RVESetHistory(Mesh,RVEtype,micro_tmp,comsolPort)

% Set correct history variable as applicable, according to RVE type
        switch RVEtype
            
            case 'ELAST1'
                
            case 'GDDAM1'    
                
            case 'PFDAM1'
                % Set history to previous converged state value
                % Required if convergence is oscillatory
                parfor idx = 1:Mesh.noElements
                    t = getCurrentTask();
                    taskID=t.ID;
                    try
                        mphstart(comsolPort(taskID));
                    catch
                    end 
                    import('com.comsol.model.*');
                    import('com.comsol.model.util.*');
                    RVEModifySolPFDAM1(idx,micro_tmp{1,idx},'history');
                end
                
            otherwise
                error('Wrong Input!')
        end

end

