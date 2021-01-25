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

function Mesh = gmsh_read(input_file_name, dir, ProblemData, ElementData)

%disp('- Reading gmsh file')

% Check number of arguments
%if nargin ~= 4
%  error('gmsh_read() takes exactly four arguments')
%endif

% Open gmsh file
file = [dir '/' input_file_name];
file_input   = fopen(file,'rt');
% ------------- Test input file ----------------------

% Read 2 first lines and check if we have mesh format 2
line = fgetl(file_input);
mesh_format = sscanf(line, '%s');
line = fgetl(file_input);
format = sscanf(line, '%f');
comp = '$MeshFormat';

%if( (strcmp(mesh_format, comp) == 0) || (format ~= 2) )
if( (strcmp(mesh_format, comp) == 0) || (strcmp(format,2) == 1) )
    error('The mesh format is NOT version 2!')
end

% ------------ Nodes --------------------

% Process file until we get the line with the number of nodes
buf = sscanf(line, '%s');
comp = '$Nodes';
while(strcmp(comp,buf) == 0)
  line = fgetl(file_input);
  buf = sscanf(line, '%s');
end

% Extract number of nodes
line = fgetl(file_input);
noNodes = sscanf(line, '%d');

% Initialise nodes matrix [x1, y1, z1; x2, y2, z2; ....; xn, yn, zn]
x = zeros(noNodes,3);

% Get nodal coordinates
for (i=1:noNodes)
  line=fgetl(file_input);
  buf = sscanf(line,'%f %f %f %f');
  x(i,:) = buf(2:4)'; %we throw away the node numbers!
end

% ------------ Elements --------------------

% Process file until we get the line with the number of elements
  comp = '$Elements';

while(strcmp(comp,buf) == 0)
  line = fgetl(file_input);
  buf = sscanf(line, '%s');    
end

% Extract number of elements
line = fgetl(file_input);
noElements = sscanf(line, '%d');

% Get first line of connectivity
line = fgetl(file_input);
buf = sscanf(line, '%d');

% Number of nodes per element
no_nodes_per_elem = size(buf,1) - (3 + buf(3));

% Get type of element
type = buf(2);

% Verify that we have the correct element
% Check number of nodes
if(no_nodes_per_elem ~= ElementData.nodes)
  error('The number of nodes per element in the mesh differ from ElementData.nodes')
end

% Check element type (gmsh 2.0 manual )
if(strcmp(ElementData.type, 'Tri3'))
  if(type ~= 2)
    error('Element type is not Tri3')
  end
elseif(strcmp(ElementData.type, 'Tri6'))
  if(type ~= 9)
    error('Element type is not Tri6')
  end
elseif(strcmp(ElementData.type, 'Tet4'))
  if(type ~=4 )
    error('Element type is not Tet4')
  end
else % Default error message
  error('Element type %s is not supported', ElementData.type)
end

% Initialise connecticity matrix and write first line
connect = zeros(noElements, no_nodes_per_elem);
connect(1,:) = buf(4 + buf(3):size(buf,1))';

% Get element connectivity
% FIXME: check that the nodes on the elements are numbered correctly!
for (i=2:noElements)
  line=fgetl(file_input);
  buf = sscanf(line,'%d');
  % Only one type of elements is allowed in the mesh
  if(type ~= buf(2))
    error('More than one type of elements is present in the mesh, did you save all elements?')
  end
  connect(i,:) = buf(4 + buf(3):size(buf,1))'; % throw away element number, type and arg list
end

% ------- Clean up and close -----------------------

% Delete coordinates
if (ProblemData.SpaceDim == 1)
  x(:,2:3) = [];
elseif (ProblemData.SpaceDim == 2)
  x(:,3) = [];
end


% Close file
fclose(file_input);

% Add members to object
Mesh.noNodes = noNodes;
Mesh.noElements = noElements;
Mesh.x = x';
Mesh.connect = connect';
