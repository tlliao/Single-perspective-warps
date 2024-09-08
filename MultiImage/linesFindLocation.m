function [ idx ] = linesFindLocation( source_lines, linepoint )
% find the index of source_lines that equals to target linepoint
% source_lines: M*4   M>=1
% linepoint: 1*4 [x1,y1,x2,y2]
idx1 = source_lines(:,1)==linepoint(1);
idx2 = source_lines(:,2)==linepoint(2);
idx3 = source_lines(:,3)==linepoint(3);
idx4 = source_lines(:,4)==linepoint(4);

idx = idx1 & idx2 & idx3 & idx4;
idx = find(idx);

end

