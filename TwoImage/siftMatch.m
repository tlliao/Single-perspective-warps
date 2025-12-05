function [pts1, pts2] = siftMatch( img1, img2)
%--------------------------------------
% SIFT keypoint detection and matching.
%--------------------------------------

%fprintf('  Keypoint detection and matching...');tic;
points1 = detectSIFTFeatures(rgb2gray(img1));
points2 = detectSIFTFeatures(rgb2gray(img2));
[features1,validPoints1] = extractFeatures(rgb2gray(img1),points1);
[features2,validPoints2] = extractFeatures(rgb2gray(img2),points2);
indexPairs = matchFeatures(features1, features2, 'MatchThreshold',10);
matchedPoints1 = validPoints1(indexPairs(:,1),:);
matchedPoints2 = validPoints2(indexPairs(:,2),:);
pts1 = double(matchedPoints1.Location');
pts2 = double(matchedPoints2.Location');
%fprintf('done (%fs)\n',toc);

% delete duplicate feature match
[~,  ind1] = unique(pts1', 'rows');
[~,  ind2] = unique(pts2', 'rows');
ind = intersect(ind1, ind2);
pts1 = pts1(:, ind);
pts2 = pts2(:, ind);

end

