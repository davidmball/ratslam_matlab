function [weight] = rs_create_posecell_weights(dim, var)
%     [WEIGHT] = rs_create_posecell_weights(DIM, VAR)
%     Creates a 3D normalised distributio of size dim^3 with a 
%     variance of var.

%     Copyright (C) 2008 David Ball (d.ball@uq.edu.au) (MATLAB version)
%     Michael Milford (m.milford1@uq.edu.au) & Gordon Wyeth (g.wyeth@uq.edu.au)
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin == 1
    error('MATLAB:rs_create_posecell_weights:moreargs', 'need at least two arguments');
end

dim_centre = floor(dim/2) + 1;

% creates a 3D normal distrubtion based on the given dimension and
% variance
weight=zeros(dim, dim, dim);
for x=1:dim
    for y=1:dim
        for z=1:dim  
           weight(x,y,z) = 1/(var*sqrt(2*pi))*exp((-(x-dim_centre)^2-(y-dim_centre)^2-(z-dim_centre)^2)/(2*var^2)); 
        end
    end
end

% ensure that it is normalised
total = sum(sum(sum(weight)));
weight = weight./total;       

end