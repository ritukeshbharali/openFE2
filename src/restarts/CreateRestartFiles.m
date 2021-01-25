%-------------------------------------------------------------------------%
% CreateRestartFiles(Sol,Res,K,Mesh,RVEdata,step,time) stores the 
% macro-scale mesh, solution, residual and stiffness matrix, and the RVE 
% files for a certain 'step' in a sub-directory under '/restart'. Note that 
% only the non-zero location and values of the stiffness matrix is stored,
% i.e., in [i, j, vals] format. Later, 'spconvert' function can be used to
% import it as a sparse matrix.
%
% INPUT:  Sol            -> Macro-scale solution
%         Res            -> Macro-scale residual
%         K              -> Macro-scale stiffness matrix
%         Mesh           -> Macro-scale mesh
%         RVEdata        -> MATLAB struct containing RVE data
%         step           -> Current step number
%         time           -> Current time
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   20.01.2021
%
%-------------------------------------------------------------------------%

function CreateRestartFiles(Sol,Res,K,Mesh,RVEdata,step,time)

% Create a directory to store restart files if it does not exist
dir_restart = '.\restart';
sub_dir     = sprintf('Step%d', step);
if ~exist(fullfile(dir_restart,sub_dir), 'dir')
    mkdir(fullfile(dir_restart,sub_dir))
end

% Save solution, residual and sparse stiffness matrix
dlmwrite(fullfile(dir_restart,sub_dir,'Solution.txt'),Sol);
dlmwrite(fullfile(dir_restart,sub_dir,'Residual.txt'),Res);
[i, j, vals] = find(K);
dlmwrite(fullfile(dir_restart,sub_dir,'SparseMatrix.txt'),[i j vals], 'delimiter', '\t');

% Store mesh
dlmwrite(fullfile(dir_restart,sub_dir,'MeshCoordinates.txt'),Mesh.x, 'delimiter', '\t');
dlmwrite(fullfile(dir_restart,sub_dir,'MeshConnectivity.txt'),Mesh.connect, 'delimiter', '\t');

% Store RVE data
fileID = fopen(fullfile(dir_restart,sub_dir,'RVEdata.txt'),'w');
fprintf(fileID,'%d\n',RVEdata.SpaceDim);
fprintf(fileID,'%s\n',RVEdata.type);
fprintf(fileID,'%d\n',RVEdata.ndofs);
fclose(fileID);

% Store time and step
dlmwrite(fullfile(dir_restart,sub_dir,'Time.txt'),[step time], 'delimiter', '\t');

% Copy RVEs to the restart folder
try
copyfile( fullfile('.\output\RVEs', 'RVE*.mph'), fullfile(dir_restart,sub_dir) );
catch
    warning('Failed to store RVEs for restart!')
end

end

