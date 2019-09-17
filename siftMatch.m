function [pts1, pts2] = siftMatch( img1, img2 )
%--------------------------------------
% SIFT keypoint detection and matching.
%--------------------------------------
fprintf('  Keypoint detection and matching...');tic;
[ kp1,ds1 ] = vl_sift(single(rgb2gray(img1)),'PeakThresh', 0,'edgethresh', 500);
[ kp2,ds2 ] = vl_sift(single(rgb2gray(img2)),'PeakThresh', 0,'edgethresh', 500);
matches   = vl_ubcmatch(ds1, ds2);
fprintf('done (%fs)\n',toc);

% extract match points' position
pts1 = kp1(1:2,matches(1,:));  pts2 = kp2(1:2,matches(2,:)); 

end

