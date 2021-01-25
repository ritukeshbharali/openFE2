% Copyright (C) 2005 Garth N .Wells
%
% Compute shape functions and determinant of the Jacobian
%
%

% Format
%
% N(i)    = shape functions of node i
% dN(i,j) = derivative of shape functions i in the direction x_j
% j       = determinant of the Jacobian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Modified by: Ritukesh Bharali
%  Date: 26 September 2020
%  *Removed all elements except 'Bar1' and 'Tri3'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N, dN, j] = ShapeFunction(x, gp, type)
  

  switch(type)
    
    case 'Bar2'
      N = zeros(2,1); dN = zeros(2,1); J = zeros(1,1); 

      xi  = gp(2);

      N(1) = -xi/2 + 0.5;
      N(2) =  xi/2 + 0.5;

      dN(1,1) = -1/2; 
      dN(2,1) =  1/2;
      
    case 'Tri3' 

      N = zeros(3,1); dN = zeros(3,2);

      xi  = gp(2);
      eta = gp(3);

      N(1)   = 1.0 - xi - eta;
      N(2)   = xi;
      N(3)   = eta;

      dN(1,1) = -1.0;
      dN(1,2) = -1.0;
      dN(2,1) =  1.0;
      dN(2,2) =  0.0;
      dN(3,1) =  0.0;
      dN(3,2) =  1.0;

    case 'Tet4'
      
      N = zeros(4,1); dN = zeros(4,3);
      
      xi  = gp(2);
      eta = gp(3);
      zi  = gp(4);

      N(1)   = 1.0 - xi - eta - zi;
      N(2)   = xi;
      N(3)   = eta;
      N(4)   = zi;

      dN(1,1:3) = -1.0;
      dN(2,1) =  1.0;
      dN(3,1) =  0.0;
      dN(3,2) =  1.0;
      dN(4,3) =  1.0;

	  otherwise 

        disp('Element type not yet programmed ')

	end

  %	transform derivatives from isoparametric configuration to global config.
  
	% form Jacobian matrix
	J = x * dN;  % Note that J is the inverse of the usual J

	% determinant of Jacobain matrix
	j = det(J);
	
  if(j < 0) 
    disp('Warning: negative Jacobian')
  end

	% transform derivatives to global system
	dN = dN /J;	
