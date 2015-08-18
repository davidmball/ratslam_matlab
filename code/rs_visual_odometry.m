function [vtrans, vrot] = rs_visual_odometry(raw_image)
%     [vtrans, vrot] = rs_visual_odometry(raw_image)

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

global prev_vrot_image_x_sums;
global prev_vtrans_image_x_sums;
global IMAGE_ODO_X_RANGE;
global IMAGE_VTRANS_Y_RANGE;
global IMAGE_VROT_Y_RANGE;
global VTRANS_SCALE;
global VISUAL_ODO_SHIFT_MATCH;

FOV_DEG = 50;
dpp = FOV_DEG / size(raw_image, 2);

% vtrans 
sub_image = raw_image(IMAGE_VTRANS_Y_RANGE, IMAGE_ODO_X_RANGE);

image_x_sums = sum(sub_image);
avint = sum(image_x_sums) / size(image_x_sums, 2);
image_x_sums = image_x_sums/avint;

[minoffset, mindiff] = rs_compare_segments(image_x_sums, prev_vtrans_image_x_sums, VISUAL_ODO_SHIFT_MATCH, size(image_x_sums, 2));

vtrans = mindiff * VTRANS_SCALE;

% a hack to detect excessively large vtrans
if vtrans > 10
    vtrans = 0;
end

prev_vtrans_image_x_sums = image_x_sums;

% now do rotation
sub_image = raw_image(IMAGE_VROT_Y_RANGE, IMAGE_ODO_X_RANGE);

image_x_sums = sum(sub_image);
avint = sum(image_x_sums) / size(image_x_sums, 2);
image_x_sums = image_x_sums/avint;

[minoffset, mindiff] = rs_compare_segments(image_x_sums, prev_vrot_image_x_sums, VISUAL_ODO_SHIFT_MATCH, size(image_x_sums, 2)); %#ok<NASGU>

vrot = minoffset * dpp * pi / 180;

prev_vrot_image_x_sums = image_x_sums;

end






