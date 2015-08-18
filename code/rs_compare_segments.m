function [offset, sdif] = rs_compare_segments(seg1, seg2, slen, cwl)
%     [offset, sdif] = rs_compare_segments(seg1, seg2, slen, cwl)
%     determine the best match between seg1 and seg2 of size cw1 allowing 
%     for shifts of up to slen

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

% assume a large difference
mindiff = 1e6;

diffs = zeros(slen);

% for each offset sum the abs difference between the two segments
for offset = 0:slen
    cdiff = abs(seg1(1 + offset:cwl) - seg2(1:cwl - offset));
    cdiff = sum(cdiff) / (cwl - offset);
    diffs(slen - offset + 1) = cdiff;
    if (cdiff < mindiff)
        mindiff = cdiff;
        minoffset = offset;
    end
end

for offset = 1:slen
    cdiff = abs(seg1(1:cwl - offset) - seg2(1 + offset:cwl));
    cdiff = sum(cdiff) / (cwl - offset);
    diffs(slen + 1 + offset) = cdiff;
    if (cdiff < mindiff)
        mindiff = cdiff;
        minoffset = -offset;
    end
end

offset = minoffset;
sdif = mindiff;

end