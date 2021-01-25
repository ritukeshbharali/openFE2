% Copyright (C) 2007 Kristian Oelgaard <k.b.oelgaard@tudelft.nl>
%
% This function processes a gmsh file (mesh format version 2) and returns an object 'Mesh'
% that is used by FECode.m developed for the class CT5123 at TU Delft. The object 'Mesh'
% has the following members:
%
% Mesh.noNodes      - the number of nodes
% Mesh.noElements   - the number of elements
% Mesh.x            - the coordinates of the nodes
% Mesh.connect      - the element connectivity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Modified by: Ritukesh Bharali
%  Date: 11.01.2021
%  *Replaced gmsh reading functionality with COMSOL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mesh = COMSOL_read(input_file_name, dir, ProblemData, ElementData)

% Open COMSOL file
file = [dir '/' input_file_name];
file_input   = fopen(file,'rt');
% ------------- Test input file ----------------------

tf = 0;
line  = fgetl(file_input);
sline = strsplit(line);
tf = strcmp('COMSOL',sline(4));
if tf ~= 1
    error('Not a COMSOL mesh!')
end

tf = 0;
while tf < 1
    line  = fgetl(file_input);
    sline = strsplit(line);
    tf = strcmp('points',sline(end));
end
noNodes = sscanf(line, '%d');

% Initialise nodes matrix [x1, y1, z1; x2, y2, z2; ....; xn, yn, zn]
x = zeros(noNodes,ProblemData.SpaceDim);

% Get nodal coordinates
tf = 0;
while tf < 1
    line  = fgetl(file_input);
    sline = strsplit(line);
    tf = strcmp('coordinates',sline(end));
end
for i=1:noNodes
  line=fgetl(file_input);
  buf = sscanf(line,'%f %f');
  x(i,:) = buf'; %we throw away the node numbers!
end

tf = 0;
while tf < 1
    line  = fgetl(file_input);
    sline = strsplit(line);
    tf = strcmp('elements',sline(end));
end
noElements = sscanf(line, '%d');

% Initialise nodes matrix [x1, y1, z1; x2, y2, z2; ....; xn, yn, zn]
connect = zeros(noElements,ElementData.nodes);

% Get elements
line  = fgetl(file_input);
for i=1:noElements
  line=fgetl(file_input);
  buf = sscanf(line,'%d %d %d');
  connect(i,:) = buf'; %we throw away the node numbers!
end

connect = connect + ones(noElements,3);

% Close file
fclose(file_input);

% Add members to object
Mesh.noNodes = noNodes;
Mesh.noElements = noElements;
Mesh.x = x';
Mesh.connect = connect';
