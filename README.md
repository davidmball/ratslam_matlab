RatSLAM Code Release (Matlab Version)

Copyright (C) 2008 David Ball (davidmichaelball@gmail.com) (MATLAB version)
Algorithm - Michael Milford & Gordon Wyeth
 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


** Introduction **

The RatSLAM system performs vision based SLAM using a computational model of the rodent hippocampus. RatSLAM is capable of performing real-time on-line SLAM in indoor and outdoor environments. The MATLAB implementation is intended to demonstrate how RatSLAM works and allow you to test your own datasets. Unlike our C version, the code is unoptimised and will slow down linerly as more visual templates and expereiences are learnt.

We are very interested in hearing about your experiences using the RatSLAM system.

The core of the RatSLAM system is:
1. Visual Template - Use normalised grayscale intensity matching to generate view templates
2. Visual Odometry or pre recorded Odometry - Use the robots wheel encoders or the change in the intensity field to get raw odometric information.
3. Pose cell attractor network - The odometric infomration shifts the activity packet(s) and view templates inject energy into the pose cells.
4. Experience Map - Topological map of expereiences with odometric links

If you use this code please reference the following paper:
David Ball, Scott Heath, Michael Milford, Gordon Wyeth, Janet Wiles (2010) A navigating rat animat, Proceedings of the 12th International Conference on the Synthesis and Simulation of Living Systems, 804-811, MIT Press

For further reading:
http://www.davidmichaelball.com/portfolio-items/openratslam/

M. J. Milford, G. Wyeth, "Mapping a Suburb with a Single Camera using a Biologically Inspired SLAM System", accepted to IEEE Transactions on Robotics Special Issue on Visual SLAM. Scheduled for publication October 2008. 

M. J. Milford, "Robot Navigation from Nature", Springer-Verlag, March 2008.

M. J. Milford, G. Wyeth, D. Prasser, "RatSLAM on the Edge: Revealing a Coherent Representation from an Overloaded Rat Brain", International Conference on Intelligent Robots and Systems, Beijing, China, 2006.

Contact davidmichaelball@gmail.com regarding the MATLAB code.


** Installation and Setup **

Either change the MATLAB path to the code directory or add it to MATLAB's path list.

Download the desired dataset(s) from https://wiki.qut.edu.au/display/cyphy/RatSLAM+MATLAB. For the xvid compressed video you may need a codec. For Windows K-Lite Codec (http://www.free-codecs.com/download/K_lite_codec_pack.htm) is good, otherwise try the ffdshow (http://sourceforge.net/projects/ffdshow-tryout/) from sourceforge.

Open the st_lucia.m or axon5.m script file and set the MOV_NAME with the full path to where you put the dataset. Run the script file.

BLOCK_READ
RENDER_RATE

** Creating datasets **

1. Create a new .m file
2. For visual odometry copy the st_lucia.m file
 and for recorded odometry copy the axon5.m file
3. Modify MOV_NAME, START_FRAME, and END_FRAME to appropriate values
4. Modify IMAGE_X_OFFSET to appropriate values
5. Specify the window for the visual templates
7. Adjust view template learning rate


** Playback **

Playing back the Experience Map, active view templates and experience templates.
The rs_exp_playback.m file allows you to playback experiences if you logged them. They will be logged at the BLOCK_READ rate in rs_main.m

Example usage is:
rs_exp_playback('st_lucia1_log', 1000 : 100 : 5000)
which will read the log files that start with st_lucia1_log starting at frame 1000 and ending at frame 5000 in increments of 100 frames. The smallest increament you can have is in BLOCK_READ's which defaults to 100.


** Advanced Parameters **

There are parameters in the rs_main.m file that will affect the core RatSLAM system. They are commented in the code itself.

Take particular note of the:
* visual template viewport
* if using visual odometry its viewport
* pose cell dimensions and exhibition/inhibition parameters.
* experience map correction loops



