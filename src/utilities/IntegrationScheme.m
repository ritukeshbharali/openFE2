% Copyright (C) 2005 Garth N .Wells
%
% Return integration point locations and weights 
%
%
% Format
%
% gp(1,j)   = integration weight for point j
% gp(2,j)   = \xi location for point j
% gp(3,j)   = \eta location for point j (2D & 3D)
% gp(4,j)   = \zi location for point j (3D only)

function gp = IntegrationScheme(points, dim)

  if dim == 1

    if points == 1	% 1-point Gauss integration in 1D

      gp     = zeros(2,1);

      gp(1) = 2.0;
      gp(2) = 0.0;

    else

      error('Only single integration point allowed!')

    end 

  elseif dim == 2

    if points == 1	% 1-point Gauss integration in 2D for triangle

      gp     = zeros(3,1);

      gp(1,1) = 0.5;                          % weight
      gp(2,1) = 1.0/3.0; gp(3,1) = 1.0/3.0;   % position
      
    else

      error('Only single integration point allowed!')

    end

  elseif dim == 3

      if points == 1	% 1-point Gauss integration in 3D for tetrahedron

      gp     = zeros(4,1);

      gp(1,1) = 1/6;                          % weight
      gp(2,1) = 0.25; 
      gp(3,1) = 0.25;
      gp(4,1) = 0.25;
      
      else

      error('Only single integration point allowed!')

      end

  else

      error('Check space dimension of the problem!')

  end

  end	

