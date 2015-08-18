function [vt_id]=rs_visual_template(raw_image, x, y, th)
% 

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

global numvts;
global vt;
global vt_history;
global prev_vt_id;
global IMAGE_VT_Y_RANGE;
global IMAGE_VT_X_RANGE;
global VT_GLOBAL_DECAY;
global VT_ACTIVE_DECAY;
global VT_SHIFT_MATCH;
global VT_MATCH_THRESHOLD;

sub_image = raw_image(IMAGE_VT_Y_RANGE, IMAGE_VT_X_RANGE);

% normalised intensity sums 
image_x_sums = sum(sub_image);
image_x_sums = image_x_sums / sum(image_x_sums);

min_offset = ones(numvts, 1);
min_diff = ones(numvts, 1);

for k = 1:numvts
    vt(k).template_decay = vt(k).template_decay - VT_GLOBAL_DECAY;
    if vt(k).template_decay < 0
        vt(k).template_decay = 0;
    end
    [min_offset(k), min_diff(k)] = rs_compare_segments(image_x_sums, vt(k).template, VT_SHIFT_MATCH, size(image_x_sums, 2));
end

[diff, diff_id] = min(min_diff);

% if this intensity template doesn't match any of the existing templates
% then create a new template
if (diff * size(image_x_sums, 2) > VT_MATCH_THRESHOLD)
    numvts = numvts + 1;
    vt(numvts).id = numvts;
    vt(numvts).template = image_x_sums;
    vt(numvts).template_decay = VT_ACTIVE_DECAY;
    vt(numvts).x_pc = x;
    vt(numvts).y_pc = y;
    vt(numvts).th_pc = th;
    vt(numvts).first = 1;           % don't want to inject energy as the vt is been created
    vt(numvts).numexps = 0;
    vt(numvts).exps = [];
    vt_id = numvts;
else
    vt_id = diff_id;
    vt(vt_id).template_decay = vt(vt_id).template_decay + VT_ACTIVE_DECAY;
    if prev_vt_id ~= vt_id
        vt(vt_id).first = 0;
    end
end

vt_history = [vt_history; vt_id];

end