function rs_main(MOV_FILE, ODO_FILE, LOG_FILE, varargin)
%     rs_main

%     Parameters could be passed better using the MATLAB method of passing
%     variables by string to overide defaults.

%     This could be much better written using Matlab's classes, especially
%     to encapsulate the globals. However, for this release we wanted to
%     support those that have older versions of Matlab.

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


% Pose cell activity network constants
% Note that increasing the dimensions will dramatically increate pose cell
% activity computation time
global PC_DIM_XY;           % The pose cell activity (x,y) dimension 
global PC_DIM_TH;           % The pose cell activity theta dimension
global PC_W_E_DIM;          % The local excitiation weight cube dimension
global PC_W_I_DIM;          % The local inhibition weight cube dimension
global PC_GLOBAL_INHIB;     % The global inhibition value
global PC_VT_INJECT_ENERGY; % The amount of energy injected when a view template is re-seen
PC_VT_INJECT_ENERGY = 0.1; 
PC_DIM_XY = 61;
PC_DIM_TH = 36;
PC_W_E_VAR = 1;             % The variance for the local excitiation guassian sphere
PC_W_E_DIM = 7;
PC_W_I_VAR = 2;             % The variance for the local inhibition guassian sphere
PC_W_I_DIM = 5;
PC_GLOBAL_INHIB = 0.00002;

% the posecell excitation and inhibition 3D weight matricies
global PC_W_EXCITE;
global PC_W_INHIB;
PC_W_EXCITE = rs_create_posecell_weights(PC_W_E_DIM, PC_W_E_VAR);
PC_W_INHIB = rs_create_posecell_weights(PC_W_I_DIM, PC_W_I_VAR);

% convienience constants
global PC_W_E_DIM_HALF;
global PC_W_I_DIM_HALF;
global PC_C_SIZE_TH;
PC_W_E_DIM_HALF = floor(PC_W_E_DIM/2);
PC_W_I_DIM_HALF = floor(PC_W_I_DIM/2);
PC_C_SIZE_TH = (2*pi)/PC_DIM_TH;

% these act as lookups to wrap the pose cell excitation/inhibition weight steps
global PC_E_XY_WRAP;
global PC_E_TH_WRAP;
global PC_I_XY_WRAP;
global PC_I_TH_WRAP;
PC_E_XY_WRAP = [(PC_DIM_XY - PC_W_E_DIM_HALF + 1):PC_DIM_XY 1:PC_DIM_XY 1:PC_W_E_DIM_HALF];
PC_E_TH_WRAP = [(PC_DIM_TH - PC_W_E_DIM_HALF + 1):PC_DIM_TH 1:PC_DIM_TH 1:PC_W_E_DIM_HALF];
PC_I_XY_WRAP = [(PC_DIM_XY - PC_W_I_DIM_HALF + 1):PC_DIM_XY 1:PC_DIM_XY 1:PC_W_I_DIM_HALF];
PC_I_TH_WRAP = [(PC_DIM_TH - PC_W_I_DIM_HALF + 1):PC_DIM_TH 1:PC_DIM_TH 1:PC_W_I_DIM_HALF];            

% these are the lookups for finding the centre of the posecell Posecells by
% rs_get_posecell_xyth()
global PC_XY_SUM_SIN_LOOKUP;
global PC_XY_SUM_COS_LOOKUP;
global PC_TH_SUM_SIN_LOOKUP;
global PC_TH_SUM_COS_LOOKUP;
global PC_CELLS_TO_AVG;
global PC_AVG_XY_WRAP;
global PC_AVG_TH_WRAP;
PC_XY_SUM_SIN_LOOKUP = sin((1:PC_DIM_XY).*2*pi/PC_DIM_XY);
PC_XY_SUM_COS_LOOKUP = cos((1:PC_DIM_XY).*2*pi/PC_DIM_XY);
PC_TH_SUM_SIN_LOOKUP = sin((1:PC_DIM_TH).*2*pi/PC_DIM_TH);
PC_TH_SUM_COS_LOOKUP = cos((1:PC_DIM_TH).*2*pi/PC_DIM_TH);
PC_CELLS_TO_AVG = 3;
PC_AVG_XY_WRAP = [(PC_DIM_XY - PC_CELLS_TO_AVG + 1):PC_DIM_XY 1:PC_DIM_XY 1:PC_CELLS_TO_AVG];
PC_AVG_TH_WRAP = [(PC_DIM_TH - PC_CELLS_TO_AVG + 1):PC_DIM_TH 1:PC_DIM_TH 1:PC_CELLS_TO_AVG];

% set the initial position in the pose network
x_pc = floor(PC_DIM_XY/2)+1; 
y_pc = floor(PC_DIM_XY/2)+1; 
th_pc = floor(PC_DIM_TH/2)+1;
global Posecells;
Posecells = zeros(PC_DIM_XY, PC_DIM_XY, PC_DIM_TH);
Posecells(x_pc, y_pc, th_pc) = 1;
global max_act_xyth_path;
max_act_xyth_path = [x_pc y_pc th_pc];

% set the initial position in the odo and experience map
odo = [0 0 pi / 2];

% specify the movie and the frames to read
movinfo = VideoReader(MOV_FILE);
START_FRAME = 1;
END_FRAME = movinfo.NumberOfFrames;

% these are the raw image dimensions
% the offset is the number of pixels from the centre of the image to the
% true zero rotation direction
IMAGE_Y_SIZE = movinfo.Height;
IMAGE_X_SIZE = movinfo.Width;

% set up the visual template module
global numvts;
global vt;
global vt_history;
global prev_vt_id;

global IMAGE_VT_Y_RANGE;
global IMAGE_VT_X_RANGE;
global VT_GLOBAL_DECAY;     % This is subtracted for all the view templates at each time step
global VT_ACTIVE_DECAY;     % This is added the best matching view template
global VT_SHIFT_MATCH;      % This determines how many +- pixels (therefore rotation) will be tested for match
global VT_MATCH_THRESHOLD;  % This threshold determines whether a new view template is generated

% The first two parameters affect the fall off of view template injection
% into the pose cell activity network
VT_GLOBAL_DECAY = 0.1;          
VT_ACTIVE_DECAY = 1.0; 
VT_SHIFT_MATCH = 20;
VT_MATCH_THRESHOLD = 0.09;
IMAGE_VT_Y_RANGE = 1:IMAGE_Y_SIZE;
IMAGE_VT_X_RANGE = 1:IMAGE_X_SIZE;

numvts = 1;
prev_vt_id = 1;
vt(1).id = 1;
vt(numvts).template_decay = 1.0;
vt(1).template(1) = 1;
vt(1).x_pc = x_pc;
vt(1).y_pc = y_pc;
vt(1).th_pc = th_pc;
vt(1).first = 1;
vt(1).numexps = 1;
vt(1).exps(1).id = 1;

vt_history = [0]; %#ok<NBRAK>

% set up the visual odometry 
global IMAGE_ODO_X_RANGE;
global IMAGE_VTRANS_Y_RANGE;
global IMAGE_VROT_Y_RANGE;
global VTRANS_SCALE;
global VISUAL_ODO_SHIFT_MATCH;

VTRANS_SCALE = 100;
VISUAL_ODO_SHIFT_MATCH = 140;
IMAGE_VTRANS_Y_RANGE = 1:IMAGE_Y_SIZE;
IMAGE_VROT_Y_RANGE = 1:IMAGE_Y_SIZE;
IMAGE_ODO_X_RANGE = 1:IMAGE_X_SIZE;

global prev_vrot_image_x_sums;
global prev_vtrans_image_x_sums;


global accum_delta_x;
global accum_delta_y;
global accum_delta_facing;
accum_delta_x = 0;
accum_delta_y = 0;
accum_delta_facing = pi/2;

global numexps;
global curr_exp_id;
global exps;
global exp_history;
global EXP_CORRECTION;          % The amount to correct each experience on either side of a link ( >0.5 is unstable)
global EXP_LOOPS;               % The number of times to run the experience map correction per frame
global EXP_DELTA_PC_THRESHOLD;  % The threshold change in pose cell activity to generate a new exp given the same view template
EXP_DELTA_PC_THRESHOLD = 1.0;
EXP_CORRECTION = 0.5;
EXP_LOOPS = 100;
exp_history = 1;

numexps = 1;
curr_exp_id = 1;
exps(1).x_pc = x_pc;
exps(1).y_pc = y_pc;
exps(1).th_pc = th_pc;
exps(1).x_m = 0;
exps(1).y_m = 0;
exps(1).facing_rad = pi/2;
exps(1).vt_id = 1;
exps(1).numlinks = 0;
exps(1).links = [];

ODO_ROT_SCALING = pi/180/7;
POSECELL_VTRANS_SCALING = 1;

% Process the parameters
for i=1:(nargin-3)
    if ischar(varargin{i})
        switch varargin{i}
            case 'RENDER_RATE', RENDER_RATE = varargin{i+1};
            case 'BLOCK_READ', BLOCK_READ = varargin{i+1};
            case 'START_FRAME', START_FRAME = varargin{i+1};
            case 'END_FRAME', END_FRAME = varargin{i+1};
                
            case 'PC_VT_INJECT_ENERGY', PC_VT_INJECT_ENERGY = varargin{i+1};               
            case 'IMAGE_VT_Y_RANGE', IMAGE_VT_Y_RANGE = varargin{i+1};
            case 'IMAGE_VT_X_RANGE', IMAGE_VT_X_RANGE = varargin{i+1};       
            case 'VT_SHIFT_MATCH', VT_SHIFT_MATCH = varargin{i+1};
            case 'VT_MATCH_THRESHOLD', VT_MATCH_THRESHOLD = varargin{i+1};

            case 'VTRANS_SCALE', VTRANS_SCALE = varargin{i+1};
            case 'VISUAL_ODO_SHIFT_MATCH', VISUAL_ODO_SHIFT_MATCH = varargin{i+1};
            case 'IMAGE_VTRANS_Y_RANGE', IMAGE_VTRANS_Y_RANGE = varargin{i+1};
            case 'IMAGE_VROT_Y_RANGE', IMAGE_VROT_Y_RANGE = varargin{i+1};
            case 'IMAGE_ODO_X_RANGE', IMAGE_ODO_X_RANGE = varargin{i+1};         
  
            case 'EXP_DELTA_PC_THRESHOLD', EXP_DELTA_PC_THRESHOLD = varargin{i+1};
            case 'EXP_CORRECTION', EXP_CORRECTION = varargin{i+1};
            case 'EXP_LOOPS', EXP_LOOPS = varargin{i+1};
                
            case 'ODO_ROT_SCALING', ODO_ROT_SCALING = varargin{i+1};
            case 'POSECELL_VTRANS_SCALING', POSECELL_VTRANS_SCALING = varargin{i+1};
        end
    end
end

vt(1).template = zeros(1, size(IMAGE_VT_X_RANGE, 2));
prev_vrot_image_x_sums = zeros(1, size(IMAGE_ODO_X_RANGE, 2));
prev_vtrans_image_x_sums = zeros(1, size(IMAGE_ODO_X_RANGE, 2));

% grab the video info and first block ... send to the vision module
mov = VideoReader(MOV_FILE); % SA replacement for aviread 234->243
vidWidth = mov.Width;
vidHeight = mov.Height;
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);
k = START_FRAME;
while  k<=(START_FRAME + (BLOCK_READ-1))
    s(k).cdata = readFrame(mov);
    k = k+1;
end

if ODO_FILE ~= 0
    ododata = csvread(ODO_FILE, START_FRAME, 0, [START_FRAME 0 min([(BLOCK_READ - 1)+START_FRAME, movinfo.NumberOfFrames]) 1]);
end

time_delta_s = [];
tic

for frame=2:min([END_FRAME, movinfo.NumberOfFrames])

    % save the experience map information to the disk for later playback
    % read the avi file in blocks and record the delta time
    if (mod(frame, BLOCK_READ) == 0)
        save(strcat(LOG_FILE, num2str(frame)), 'frame', 'exps', 'exp_history', 'vt_history');
        time_delta_s = [time_delta_s; toc]; %#ok<AGROW>
        while  k<= min([(frame+BLOCK_READ - 1)+START_FRAME, movinfo.NumberOfFrames]);
         s(k).cdata = readFrame(mov);
         k = k+1;
         end
        if ODO_FILE ~= 0
            ododata = csvread(ODO_FILE, frame+START_FRAME, 0, [frame+START_FRAME 0 min([(frame+BLOCK_READ - 1)+START_FRAME, movinfo.NumFrames]) 1]);
        end
        tic
    end

    % visual templates and visual odo uses intensity so convert to grayscale
    im = rgb2gray(s(mod(frame, BLOCK_READ) + 1).cdata);

    % get the most active view template
    [vt_id] = rs_visual_template(im, x_pc, y_pc, th_pc);
    
    % get the odometry from the video
    if ODO_FILE == 0
        [vtrans, vrot] = rs_visual_odometry(im);
    else
        vtrans = ododata(mod(frame, BLOCK_READ) + 1, 1);
        vrot = ododata(mod(frame, BLOCK_READ) + 1, 2)*ODO_ROT_SCALING;
    end
    % track the motion based on raw odometry only (for comparison)
    odo(frame, 1) = odo(frame - 1, 1) + vtrans * cos(odo(frame - 1, 3) + vrot);
    odo(frame, 2) = odo(frame - 1, 2) + vtrans * sin(odo(frame - 1, 3) + vrot);
    odo(frame, 3) = odo(frame - 1, 3) + vrot;

    % update the pose cells accounts for the active vt and PI
    rs_posecell_iteration(vt_id, vtrans * POSECELL_VTRANS_SCALING, vrot);

    % gets the 'best' centre of the pose cell activity
    [x_pc, y_pc, th_pc] = rs_get_posecell_xyth();
 
    % runs an interation of the experience map
    rs_experience_map_iteration(vt_id, vtrans, vrot, x_pc, y_pc, th_pc);
    
    % render debug information
    if (mod(frame, RENDER_RATE) == 0)
        
        % render the raw image
        subplot(3, 3, 1, 'replace');
        image(s(mod(frame, BLOCK_READ) + 1).cdata);
        title('Raw Image');
        
        % render the history of visual templates
        subplot(3, 3, 2, 'replace');
        plot(vt_history, '.');
        title('Frame vs View Template');

        % render the agent path given by raw odometry
        subplot(3, 3, 3, 'replace');
        plot(odo(:, 1), odo(:, 2), 'x');
        axis equal;
        title('Raw Odometry');
        
%         % render the time to process each block
%         subplot(3, 3, 5, 'replace');
%         plot(time_delta_s);
%         title('Time versus CPU time');
        
        % render the pose cell network including the xy path
        subplot(3, 3, 4, 'replace');
        max_act_xyth_path = [max_act_xyth_path; x_pc y_pc th_pc]; %#ok<AGROW>
        phandles = contourslice(Posecells, 1:4:PC_DIM_XY, 1:4:PC_DIM_XY, 1:3:PC_DIM_TH, 3);
        axis([1 PC_DIM_XY 1 PC_DIM_XY 1 PC_DIM_TH]);
        view(3);
        set(phandles,'LineWidth',0.1);
        grid on;
        hold on;
        plot3(max_act_xyth_path(:,2), max_act_xyth_path(:,1), 0 * max_act_xyth_path(:,3), '.m');
        plot3([max_act_xyth_path(end,2) max_act_xyth_path(end,2)], [max_act_xyth_path(end,1) max_act_xyth_path(end,1)], [0 th_pc]);
        hold off;
        title('Pose Cell Activity');
        
        % render the history of experiences
        subplot(3, 3, 7, 'replace');
        plot(exp_history, '.');
        title('Frame vs Experience');

        % render the experience map
        subplot(2, 2, 4, 'replace');
        plot([exps(:).x_m], [exps(:).y_m], 'x');
        axis equal;
        title('Experience Map');
        
        drawnow;

    end
    
    prev_vt_id = vt_id;
end

total_time = sum(time_delta_s) %#ok<NOPRT,NASGU,NOPTS>

end