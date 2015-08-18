function [X, Y, TH] = rs_get_posecell_xyth()
%     [X, Y, TH] = rs_get_posecell_xyth()
%     Returns the approximate averaged centre of the most active activity
%     packet. This implementation averages the cells around the maximally
%     activated cell.

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

global PC_DIM_XY;
global PC_DIM_TH;
global PC_XY_SUM_SIN_LOOKUP;
global PC_XY_SUM_COS_LOOKUP;
global PC_TH_SUM_SIN_LOOKUP;
global PC_TH_SUM_COS_LOOKUP;
global PC_CELLS_TO_AVG;
global PC_AVG_XY_WRAP;
global PC_AVG_TH_WRAP;

global Posecells;


% find the max activated cell
indexes = find(Posecells);
[value, index] = max(Posecells(indexes));
[x y z] = ind2sub(size(Posecells), indexes(index));

% take the max activated cell +- AVG_CELL in 3d space
z_Posecells=zeros(PC_DIM_XY, PC_DIM_XY, PC_DIM_TH);
z_Posecells(PC_AVG_XY_WRAP(x:x+PC_CELLS_TO_AVG*2), PC_AVG_XY_WRAP(y:y+PC_CELLS_TO_AVG*2), PC_AVG_TH_WRAP(z:z+PC_CELLS_TO_AVG*2)) = ...
    Posecells(PC_AVG_XY_WRAP(x:x+PC_CELLS_TO_AVG*2), PC_AVG_XY_WRAP(y:y+PC_CELLS_TO_AVG*2), PC_AVG_TH_WRAP(z:z+PC_CELLS_TO_AVG*2));

% get the sums for each axis
x_sums = sum(sum(z_Posecells, 2), 3)';
y_sums = sum(sum(z_Posecells, 1), 3);
th_sums = sum(sum(z_Posecells, 1), 2);
th_sums = th_sums(:)';

% now find the (x, y, th) using population vector decoding to handle the wrap around 
X = mod(atan2(sum(PC_XY_SUM_SIN_LOOKUP.*x_sums), sum(PC_XY_SUM_COS_LOOKUP.*x_sums))*PC_DIM_XY/(2*pi), PC_DIM_XY);
Y = mod(atan2(sum(PC_XY_SUM_SIN_LOOKUP.*y_sums), sum(PC_XY_SUM_COS_LOOKUP.*y_sums))*PC_DIM_XY/(2*pi), PC_DIM_XY);
TH = mod(atan2(sum(PC_TH_SUM_SIN_LOOKUP.*th_sums), sum(PC_TH_SUM_COS_LOOKUP.*th_sums))*PC_DIM_TH/(2*pi), PC_DIM_TH);

end