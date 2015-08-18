function rs_experience_map_iteration(vt_id, vtrans, vrot, x_pc, y_pc, th_pc)
%     rs_experience_map_iteration

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

global exps;
global vt;
global prev_vt_id;
global curr_exp_id;
global numexps;
global prev_exp_id;
global accum_delta_x;
global accum_delta_y;
global accum_delta_facing;

global PC_DIM_XY;
global PC_DIM_TH;

global EXP_DELTA_PC_THRESHOLD

% integrate the delta x, y, facing
accum_delta_facing = Clip_Rad_180(accum_delta_facing + vrot);
accum_delta_x = accum_delta_x + vtrans * cos(accum_delta_facing);
accum_delta_y = accum_delta_y + vtrans * sin(accum_delta_facing);

delta_pc = sqrt(Get_Min_Delta(exps(curr_exp_id).x_pc, x_pc, PC_DIM_XY)^2 + Get_Min_Delta(exps(curr_exp_id).y_pc, y_pc, PC_DIM_XY)^2 + Get_Min_Delta(exps(curr_exp_id).th_pc, th_pc, PC_DIM_TH)^2);

% if the vt is new or the pc x,y,th has changed enough create a new
% experience
if vt(vt_id).numexps == 0 || delta_pc > EXP_DELTA_PC_THRESHOLD

    numexps = numexps + 1;
    Create_New_Exp(curr_exp_id, numexps, x_pc, y_pc, th_pc, vt_id);

    prev_exp_id = curr_exp_id;
    curr_exp_id = numexps;

    accum_delta_x = 0;
    accum_delta_y = 0;
    accum_delta_facing = exps(curr_exp_id).facing_rad;

    % if the vt has changed (but isn't new) search for the matching exp
elseif vt_id ~= prev_vt_id

    % find the exp associated with the current vt and that is under the
    % threshold distance to the centre of pose cell activity
    % if multiple exps are under the threshold then don't match (to reduce
    % hash collisions)
    matched_exp_id = 0;
    matched_exp_count = 0;
%    [x_pc y_pc th_pc vt_id]
    
    for search_id=1:vt(vt_id).numexps
%        [search_id exps(vt(vt_id).exps(search_id).id).x_pc exps(vt(vt_id).exps(search_id).id).y_pc exps(vt(vt_id).exps(search_id).id).th_pc]
        delta_pc(search_id) = sqrt(Get_Min_Delta(exps(vt(vt_id).exps(search_id).id).x_pc, x_pc, PC_DIM_XY)^2 + Get_Min_Delta(exps(vt(vt_id).exps(search_id).id).y_pc, y_pc, PC_DIM_XY)^2 + Get_Min_Delta(exps(vt(vt_id).exps(search_id).id).th_pc, th_pc, PC_DIM_TH)^2);
        if delta_pc(search_id) < EXP_DELTA_PC_THRESHOLD
           matched_exp_count = matched_exp_count + 1; 
        end
    end
        
    if matched_exp_count > 1
        % this means we aren't sure about which experience is a match due
        % to hash table collision
        % instead of a false posivitive which may create blunder links in
        % the experience map keep the previous experience
%        matched_exp_count
        
    else
        [min_delta, min_delta_id] = min(delta_pc);

        if min_delta < EXP_DELTA_PC_THRESHOLD

            matched_exp_id = vt(vt_id).exps(min_delta_id).id;

            % see if the prev exp already has a link to the current exp
            link_exists = 0;
            for link_id = 1:exps(curr_exp_id).numlinks
                if exps(curr_exp_id).links(link_id).exp_id == matched_exp_id
                    link_exists = 1;
                    break;
                end
            end

%            [matched_exp_id search_id link_exists]

            % if we didn't find a link then create the link between current
            % experience and the expereince for the current vt
            if link_exists == 0
                exps(curr_exp_id).numlinks = exps(curr_exp_id).numlinks + 1;
                exps(curr_exp_id).links(exps(curr_exp_id).numlinks).exp_id = matched_exp_id;
                exps(curr_exp_id).links(exps(curr_exp_id).numlinks).d = sqrt(accum_delta_x^2 + accum_delta_y^2);
                exps(curr_exp_id).links(exps(curr_exp_id).numlinks).heading_rad = Get_Signed_Delta_Rad(exps(curr_exp_id).facing_rad, atan2(accum_delta_y, accum_delta_x));
                exps(curr_exp_id).links(exps(curr_exp_id).numlinks).facing_rad = Get_Signed_Delta_Rad(exps(curr_exp_id).facing_rad, accum_delta_facing);
            end

        end

        % if there wasn't an experience with the current vt and posecell x y th
        % then create a new experience
        if matched_exp_id == 0
            numexps = numexps + 1;
%            numexps
            Create_New_Exp(curr_exp_id, numexps, x_pc, y_pc, th_pc, vt_id);
            matched_exp_id = numexps;
        end

        prev_exp_id = curr_exp_id;
        curr_exp_id = matched_exp_id;

        accum_delta_x = 0;
        accum_delta_y = 0;
        accum_delta_facing = exps(curr_exp_id).facing_rad;
    end
    
end

global EXP_CORRECTION;
global EXP_LOOPS;

% do the experience map correction interatively for all the links in all
% the experiences
for i=1:EXP_LOOPS

    for exp_id=1:numexps

        for link_id=1:exps(exp_id).numlinks

            % experience 0 has a link to experience 1
            e0 = exp_id;
            e1 = exps(exp_id).links(link_id).exp_id;

            % work out where e0 thinks e1 (x,y) should be based on the stored
            % link information
            lx = exps(e0).x_m + exps(e0).links(link_id).d * cos(exps(e0).facing_rad + exps(e0).links(link_id).heading_rad);
            ly = exps(e0).y_m + exps(e0).links(link_id).d * sin(exps(e0).facing_rad + exps(e0).links(link_id).heading_rad);

            % correct e0 and e1 (x,y) by equal but opposite amounts
            % a 0.5 correction parameter means that e0 and e1 will be fully
            % corrected based on e0's link information
            exps(e0).x_m = exps(e0).x_m + (exps(e1).x_m - lx) * EXP_CORRECTION;
            exps(e0).y_m = exps(e0).y_m + (exps(e1).y_m - ly) * EXP_CORRECTION;
            exps(e1).x_m = exps(e1).x_m - (exps(e1).x_m - lx) * EXP_CORRECTION;
            exps(e1).y_m = exps(e1).y_m - (exps(e1).y_m - ly) * EXP_CORRECTION;

            % determine the angle between where e0 thinks e1's facing
            % should be based on the link information
            df = Get_Signed_Delta_Rad((exps(e0).facing_rad + exps(e0).links(link_id).facing_rad), exps(e1).facing_rad);
            
            % correct e0 and e1 facing by equal but opposite amounts
            % a 0.5 correction parameter means that e0 and e1 will be fully
            % corrected based on e0's link information           
            exps(e0).facing_rad = Clip_Rad_180(exps(e0).facing_rad + df * EXP_CORRECTION);
            exps(e1).facing_rad = Clip_Rad_180(exps(e1).facing_rad - df * EXP_CORRECTION);
        end
    end

end

% keep a frame by frame history of which experience was currently active
global exp_history;
exp_history = [exp_history; curr_exp_id];

end

function Create_New_Exp(curr_exp_id, new_exp_id, x_pc, y_pc, th_pc, vt_id)
% create a new experience and current experience to it
global vt;
global exps;
global accum_delta_x;
global accum_delta_y;
global accum_delta_facing;

% add link information to the current experience for the new experience
% including the experience_id, odo distance to the experience, odo heading
% (relative to the current experience's facing) to the experience, odo delta 
% facing (relative to the current expereience's facing).
exps(curr_exp_id).numlinks = exps(curr_exp_id).numlinks + 1;
exps(curr_exp_id).links(exps(curr_exp_id).numlinks).exp_id = new_exp_id;
exps(curr_exp_id).links(exps(curr_exp_id).numlinks).d = sqrt(accum_delta_x^2 + accum_delta_y^2);
exps(curr_exp_id).links(exps(curr_exp_id).numlinks).heading_rad = Get_Signed_Delta_Rad(exps(curr_exp_id).facing_rad, atan2(accum_delta_y, accum_delta_x));
exps(curr_exp_id).links(exps(curr_exp_id).numlinks).facing_rad = Get_Signed_Delta_Rad(exps(curr_exp_id).facing_rad, accum_delta_facing);

% create the new experience which will have no links to being with
exps(new_exp_id).x_pc = x_pc;
exps(new_exp_id).y_pc = y_pc;
exps(new_exp_id).th_pc = th_pc;
exps(new_exp_id).vt_id = vt_id;
exps(new_exp_id).x_m = exps(curr_exp_id).x_m + accum_delta_x;
exps(new_exp_id).y_m = exps(curr_exp_id).y_m + accum_delta_y;
exps(new_exp_id).facing_rad = Clip_Rad_180(accum_delta_facing);
exps(new_exp_id).numlinks = 0;
exps(new_exp_id).links = [];

% add this experience id to the vt for efficient lookup
vt(vt_id).numexps = vt(vt_id).numexps + 1;
vt(vt_id).exps(vt(vt_id).numexps).id = new_exp_id;
end

function [angle]=Clip_Rad_360(angle)
% Clip the input angle to between 0 and 2pi radians

while angle < 0
    angle = angle + 2*pi;
end
while angle >= 2*pi
    angle = angle - 2*pi;
end
end

function [angle]=Clip_Rad_180(angle)
% Clip the input angle to between -pi and pi radians

while angle > pi
    angle = angle - 2*pi;
end
while angle <= -pi
    angle = angle + 2*pi;
end
end

function [delta]=Get_Min_Delta(d1, d2, max)
% Get the minimum delta distance between two values assuming a wrap to zero
% at max

delta = min([abs(d1 - d2), max - abs(d1 - d2)]);
end


function [angle]=Get_Signed_Delta_Rad(angle1, angle2)
% Get the signed delta angle from angle1 to angle2 handling the wrap from 2pi
% to 0.

dir = Clip_Rad_180(angle2 - angle1);

delta_angle = abs(Clip_Rad_360(angle1) - Clip_Rad_360(angle2));

if (delta_angle) < (2*pi - delta_angle)
    if (dir > 0)
        angle = delta_angle;
    else
        angle = -delta_angle;
    end
else
    if (dir > 0)
        angle = 2*pi - delta_angle;
    else
        angle = -(2*pi - delta_angle);
    end
end

end

