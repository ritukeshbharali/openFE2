%-------------------------------------------------------------------------%
% RVERevert(Mesh,RVEtype,micro_tmp,comsolPort) is a master function that
% calls the apprpriate function to revert the RVE solution to an earlier
% converged state 'micro_tmp'.
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

function RVERevert(Mesh,RVEtype,micro_tmp,comsolPort)

% Set correct history variable as applicable, according to RVE type
        switch RVEtype
            
            case 'ELAST1'
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
                    RVEModifySolELAST1(idx,micro_tmp{1,idx},'revert');
                end
                
            case 'GDDAM1'
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
                    RVEModifySolGDDAM1(idx,micro_tmp{1,idx},'revert');
                end                
                
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
                    RVEModifySolPFDAM1(idx,micro_tmp{1,idx},'revert');
                end
                
            otherwise
                error('Wrong Input!')
        end

end

