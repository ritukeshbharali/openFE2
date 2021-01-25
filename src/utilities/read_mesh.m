% Copyright (C) 2007 Kristian Oelgaard <k.b.oelgaard@tudelft.nl>
%
% This function processes either a a gmsh file (mesh format version 2) or a hand made mesh file
% and returns an object 'Mesh' that is used by FECode.m developed for the class CT5123 at TU Delft.
% The object 'Mesh' has the following members:
%
% Mesh.noNodes      - the number of nodes
% Mesh.noElements   - the number of elements
% Mesh.x            - the coordinates of the nodes
% Mesh.connect      - the element connectivity
%
% This function assumes that the name of the file is either 'mesh.msh' (gmsh) or 'mesh.m'.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Modified by: Ritukesh Bharali
%  Date: 26 September 2020
%  *Removed limitation on the mesh file name
%  *Calls gmsh_read or COMSOL_read based on file extension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mesh = read_mesh(dir, ProblemData)

ElementData.nodes = ProblemData.nodes;
ElementData.type  = ProblemData.type;

name_mesh = ProblemData.mesh_file;
[~,~,ext] = fileparts(name_mesh);

switch ext
    case '.msh'
        %disp('- Reading GMSH mesh file')
        Mesh = gmsh_read(name_mesh, dir, ProblemData, ElementData);
        Mesh.type = 'GMSH';
    case '.mphtxt'
        %disp('- Reading COMSOL mesh file')
        Mesh = COMSOL_read(name_mesh, dir, ProblemData, ElementData);
        Mesh.type = 'COMSOL';
    case '.m'
        MATLAB_mesh
    otherwise
        error('Mesh file not recognised!')
end

end
