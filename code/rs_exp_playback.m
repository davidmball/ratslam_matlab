function rs_exp_playback(FILENAME, FRAMES)
%     rs_exp_playback(FILENAME, FRAMES)
%     Plays back the logged experience map, visual template and experience
%     templates given by FILENAME for FRAMES (start:increment:end).

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

for frame=FRAMES

    load(strcat(FILENAME, num2str(frame)));

    subplot(2,2,1, 'replace');
    plot(vt_history, '.');
    subplot(2,2,2, 'replace');
    plot(exp_history, '.');
    subplot(2,1,2, 'replace');
    plot([exps(:).x_m], [exps(:).y_m], 'x');        
    axis equal;

    drawnow;
end

end