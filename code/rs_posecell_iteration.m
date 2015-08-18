function rs_posecell_iteration(vt_id, vtrans, vrot)
% rs_posecell_iteration

% Pose cell update steps
% 1. Add view template energy
% 2. Local excitation
% 3. Local inhibition
% 4. Global inhibition
% 5. Normalisation
% 6. Path Integration (vtrans then vrot)

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
global PC_W_E_DIM;
global PC_W_I_DIM;
global PC_GLOBAL_INHIB;
global PC_C_SIZE_TH;
global PC_E_XY_WRAP;
global PC_E_TH_WRAP;
global PC_I_XY_WRAP;
global PC_I_TH_WRAP;
global PC_W_EXCITE;
global PC_W_INHIB;
global PC_VT_INJECT_ENERGY;

global Posecells;
global vt;

% if this isn't a new vt then add the energy at its associated posecell
% location
if vt(vt_id).first ~= 1
    act_x = min([max([round(vt(vt_id).x_pc), 1]), PC_DIM_XY]);
    act_y = min([max([round(vt(vt_id).y_pc), 1]), PC_DIM_XY]);
    act_th = min([max([round(vt(vt_id).th_pc), 1]), PC_DIM_TH]);
    
    % this decays the amount of energy that is injected at the vt's
    % posecell location
    % this is important as the posecell Posecells will errounously snap 
    % for bad vt matches that occur over long periods (eg a bad matches that
    % occur while the agent is stationary). This means that multiple vt's
    % need to be recognised for a snap to happen
    energy = PC_VT_INJECT_ENERGY * 1/30 * (30 - exp(1.2 * vt(vt_id).template_decay));
    if energy > 0
        Posecells(act_x, act_y, act_th) = Posecells(act_x, act_y, act_th) + energy;
    end
end


% local excitation - PC_le = PC elements * PC weights
pca_new = zeros(PC_DIM_XY, PC_DIM_XY, PC_DIM_TH);
for x=1:PC_DIM_XY
    for y=1:PC_DIM_XY
        for z=1:PC_DIM_TH
            if Posecells(x,y,z) ~= 0
                pca_new(PC_E_XY_WRAP(x:x+PC_W_E_DIM-1),PC_E_XY_WRAP(y:y+PC_W_E_DIM-1),PC_E_TH_WRAP(z:z+PC_W_E_DIM-1)) = ...
                    pca_new(PC_E_XY_WRAP(x:x+PC_W_E_DIM-1),PC_E_XY_WRAP(y:y+PC_W_E_DIM-1),PC_E_TH_WRAP(z:z+PC_W_E_DIM-1)) + Posecells(x,y,z).*PC_W_EXCITE;
            end
        end
    end
end
Posecells = pca_new;

% local inhibition - PC_li = PC_le - PC_le elements * PC weights
pca_new = zeros(PC_DIM_XY, PC_DIM_XY,PC_DIM_TH);  
for x=1:PC_DIM_XY
    for y=1:PC_DIM_XY
        for z=1:PC_DIM_TH
            if Posecells(x,y,z) ~= 0
                pca_new(PC_I_XY_WRAP(x:x+PC_W_I_DIM-1),PC_I_XY_WRAP(y:y+PC_W_I_DIM-1),PC_I_TH_WRAP(z:z+PC_W_I_DIM-1)) = ...
                    pca_new(PC_I_XY_WRAP(x:x+PC_W_I_DIM-1),PC_I_XY_WRAP(y:y+PC_W_I_DIM-1),PC_I_TH_WRAP(z:z+PC_W_I_DIM-1)) + Posecells(x,y,z).*PC_W_INHIB;
            end
        end
    end
end
Posecells = Posecells - pca_new;

% local global inhibition - PC_gi = PC_li elements - inhibition
Posecells = (Posecells >= PC_GLOBAL_INHIB) .* (Posecells - PC_GLOBAL_INHIB);

% normalisation
total = sum(sum(sum(Posecells)));
Posecells = Posecells./total;


% Path Integration
% vtrans affects xy direction
% shift in each th given by the th
for dir_pc=1:PC_DIM_TH        

    % radians
    dir = (dir_pc - 1) * PC_C_SIZE_TH;
    
    % N,E,S,W are straightforward
    if dir == 0
        Posecells(:,:,dir_pc) = Posecells(:,:,dir_pc).*(1.0 - vtrans) + circshift(Posecells(:,:,dir_pc), [0 1]).*vtrans;
    elseif dir == pi/2
        Posecells(:,:,dir_pc) = Posecells(:,:,dir_pc).*(1.0 - vtrans) + circshift(Posecells(:,:,dir_pc), [1 0]).*vtrans;
    elseif dir == pi
        Posecells(:,:,dir_pc) = Posecells(:,:,dir_pc).*(1.0 - vtrans) + circshift(Posecells(:,:,dir_pc), [0 -1]).*vtrans;
    elseif dir == 3*pi/2
        Posecells(:,:,dir_pc) = Posecells(:,:,dir_pc).*(1.0 - vtrans) + circshift(Posecells(:,:,dir_pc), [-1 0]).*vtrans;        
    else
        % rotate the Posecells instead of implementing for four quadrants
        pca90 = rot90(Posecells(:,:,dir_pc), floor(dir *2/pi));
        dir90 = dir - floor(dir *2/pi) * pi/2;

        % extend the Posecells one unit in each direction (max supported at the moment)
        % work out the weight contribution to the NE cell from the SW, NW, SE cells 
        % given vtrans and the direction
        % weight_sw = v * cos(th) * v * sin(th)
        % weight_se = (1 - v * cos(th)) * v * sin(th)
        % weight_nw = (1 - v * sin(th)) * v * sin(th)
        % weight_ne = 1 - weight_sw - weight_se - weight_nw
        % think in terms of NE divided into 4 rectangles with the sides
        % given by vtrans and the angle
        pca_new=zeros(PC_DIM_XY+2);            
        pca_new(2:end-1,2:end-1) = pca90;
        weight_sw = vtrans^2 *cos(dir90) * sin(dir90);
        weight_se = vtrans*sin(dir90) - vtrans^2 *cos(dir90) * sin(dir90);
        weight_nw = vtrans*cos(dir90) - vtrans^2 *cos(dir90) * sin(dir90);
        weight_ne = 1.0 - weight_sw - weight_se - weight_nw;
  
        % circular shift and multiple by the contributing weight
        % copy those shifted elements for the wrap around       
        pca_new = pca_new.*weight_ne + circshift(pca_new, [0 1]).*weight_nw + circshift(pca_new, [1 0]).*weight_se + circshift(pca_new, [1 1]).*weight_sw;
        pca90 = pca_new(2:end-1,2:end-1);
        pca90(2:end,1) = pca90(2:end,1) + pca_new(3:end-1,end);
        pca90(1,2:end) = pca90(1,2:end) + pca_new(end,3:end-1);
        pca90(1,1) = pca90(1,1) + pca_new(end:end);

        % unrotate the pose cell xy layer
        Posecells(:,:,dir_pc) = rot90(pca90, 4 - floor(dir * 2/pi));
    end
end


% Path Integration - Theta
% Shift the pose cells +/- theta given by vrot
if vrot ~= 0
    % mod to work out the partial shift amount
    weight = mod(abs(vrot)/PC_C_SIZE_TH, 1);
    if weight == 0
        weight = 1.0;
    end
    Posecells = circshift(Posecells, [0 0 sign(vrot) * floor(abs(vrot)/PC_C_SIZE_TH)]) * (1.0 - weight) ...
            + circshift(Posecells, [0 0 sign(vrot) * ceil(abs(vrot)/PC_C_SIZE_TH)]) * (weight); 
end


end











