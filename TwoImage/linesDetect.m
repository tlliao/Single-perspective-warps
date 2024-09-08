function [ sample_lines, slope_lines ] = linesDetect( imgpath, img, C1, C2 )
% detect line segments in image using LSD: "A fast line segment detector with a
% false detection control" by R. G. von Gioi, J. Jakubowicz, J.-M. Morel, and G. Randall. Lsd
%% get the start_points and end_points of each straight line use LSD.
lines = lsd(imgpath);
[M, N, ~] = size(img);
x_dis = (N-1)/C2; y_dis = (M-1)/C1;  % width and height of each grid

len_threshold = sqrt(x_dis^2+ y_dis^2);  % length threshold of line segments, ignore all segments of length less than threshold
len_lines = sqrt( (lines(2,:)-lines(1,:)).^2 + (lines(4,:)-lines(3,:)).^2 );
lines_ = lines(:,len_lines>=len_threshold);
lines_(6,:) = len_lines(len_lines>=len_threshold);

%% plot the lines.
% figure,imshow(img);
% hold on
% for k = 1:size(lines, 2)
%     plot(lines(1:2, k), lines(3:4, k), 'LineWidth', lines(5, k) / 2, 'Color', [1, 0, 0]);
% end
% hold off
% % 
% figure,imshow(img);
% hold on
% for k = 1:size(lines_, 2)
%     plot(lines_(1:2, k), lines_(3:4, k), 'LineWidth', lines_(5, k)/2, 'Color', [1, 0, 0]); % / 
%     text((lines_(1,k)+lines_(2,k))/2,(lines_(3,k)+lines_(4,k))/2, num2str(k),'FontSize',12);
% end
% hold off

%% uniformly sample points on lines such that each mesh grid contains a sample point of the line
num_lines = max(3, 2*round(lines_(6,:)./min(x_dis, y_dis)));  % sample numbers of each line segment
sample_lines = zeros( 2*size(lines_,2), max(num_lines)+1);
slope_lines = zeros(2*size(lines_,2),1);  % slope of each lines
for k=1:size(lines_,2)
    x1 = lines_(1,k);  y1 = lines_(3,k);
    x2 = lines_(2,k);  y2 = lines_(4,k);
    if x1~=x2
        slope = (y2-y1)/(x2-x1);
        xseq = linspace(x1,x2, num_lines(k));
        yseq = (xseq-x1).*slope+y1;
        r_xseq = (xseq>=1) & (xseq<=N);
        r_yseq = (yseq>=1) & (yseq<=M);
        r_seq = r_xseq & r_yseq;
        sample_lines(2*k-1:2*k, 1:sum(r_seq)) = [xseq(r_seq); yseq(r_seq)];
        sample_lines(2*k-1:2*k, end) = sum(r_seq); 
        slope_lines(2*k-1:2*k) = slope;
    end
    if x1==x2
        xseq = x1.*ones(1, num_lines(k));
        yseq = linspace(y1,y2, num_lines(k));
        r_xseq = (xseq>=1) & (xseq<=N);
        r_yseq = (yseq>=1) & (yseq<=M);
        r_seq = r_xseq & r_yseq;
        sample_lines(2*k-1:2*k, 1:sum(r_seq)) = [xseq(r_seq); yseq(r_seq)];
        sample_lines(2*k-1:2*k, end) = sum(r_seq);    
        slope_lines(2*k-1:2*k) = inf;
    end
end

end

