%-------------------------------------------------------------------------%
% [K,Res,Sol,RVEdata,step,time] = ImportRestartFiles(dir_restart,sub_dir)
% reads all files in directory 'dir_restart/sub_dir' and yields the
% current macro-scale sparse stiffness matrix, residual, solution along
% with RVE data (SpaceDim,type,ndofs) and time, step information. These
% enables one to restart simulations.
%
% INPUT:  dir_restart    -> Restart directory
%         sub_dir        -> Directory containing restart files of a
%                           particular step
%
% OUTPUT: K              -> Macro-scale sparse stiffness matrix
%         Res            -> Macro-scale residual
%         Sol            -> Macro-scale solution
%         RVEdata        -> MATLAB struct containing RVE data (SpaceDim,
%                           type, ndofs)
%         step           -> Current step
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

function [K,Res,Sol,RVEdata,step,time] = ImportRestartFiles(dir_restart,sub_dir)

% Import sparse stiffness matrix, residual and solution
K   = dlmread(fullfile(dir_restart,sub_dir,'SparseMatrix.txt'),'\t');
K   = spconvert(K);
Res = dlmread(fullfile(dir_restart,sub_dir,'Residual.txt'));
Sol = dlmread(fullfile(dir_restart,sub_dir,'Solution.txt'));

% Import RVE data
fileID = fopen(fullfile(dir_restart,sub_dir,'RVEdata.txt'),'r');
RVEdata.SpaceDim = str2double(fgetl(fileID));
RVEdata.type     = fgetl(fileID);
RVEdata.ndofs    = str2double(fgetl(fileID));
fclose(fileID);

% Import step and time
time = dlmread(fullfile(dir_restart,sub_dir,'Time.txt'));
step = time(1);
time = time(2);

% Copy RVEs to the working directory './output/RVEs'
try
copyfile( fullfile(dir_restart,sub_dir,'RVE*.mph'), '.\output\RVEs'  );
catch
    warning('Failed to get RVEs for the restart!')
end

end

