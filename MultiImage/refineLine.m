function [ refine_lines ] = refineLine( lines, img )
% refine lines in image coordinates,
% lines: M*4 matrix M lines
% |x1,y1,x2,y2|;...
[sz1, sz2, ~] = size(img);
refine_lines = lines;
for i=1:size(lines,1)
    x1=lines(i,1); x2=lines(i,3);
    y1=lines(i,2); y2=lines(i,4);
    if x1>=1 && x1<=sz2 && x2>=1 && x2<=sz2 && y1>=1 && y1<=sz1 && y2>=1 && y2<=sz1
        continue;        
    end
    a=y2-y1; b=x1-x2; c=x2*y1-x1*y2;
    if abs(b)<=eps  % x1=x2
        x1=min(max(1.1, x1),sz2); x2=min(max(1.1, x2),sz2);
        y1=min(max(1.1, y1),sz1); y2=min(max(1.1, y2),sz1);
    else
        x1=min(max(1.1, x1),sz2); x2=min(max(1.1, x2),sz2);
        y1=-(a*x1+c)/b; y2=-(a*x2+c)/b;
    end
    if abs(a)>eps
        y1=min(max(1.1, y1),sz1); y2=min(max(1.1, y2),sz1);
        x1=-(b*y1+c)/a;  x2=-(b*y2+c)/a;
    end
    refine_lines(i,1)=x1; refine_lines(i,3)=x2;
    refine_lines(i,2)=y1; refine_lines(i,4)=y2;
end


end

