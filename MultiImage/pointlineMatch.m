function [data_orig, inliers, data_lines, h, imgs_pairs] = pointlineMatch( cell_imgs, cell_paths)
%--------------------------------------
% SIFT keypoint detection and matching.
%--------------------------------------
% global model specific function handlers.
clear global;
global fitfn resfn degenfn psize numpar
fitfn = 'homography_fit';
resfn = 'homography_res';
degenfn = 'homography_degen';
psize   = 4;
numpar  = 9;


img_n = size(cell_imgs,1);  % number of images
imgs_pairs = [(1:img_n-1)', (2:img_n)']; % generating image pairs with all available image idexes.
imgs_pairs = [imgs_pairs ones(size(imgs_pairs,1),1)]; % last column indicates if the pair is valid or not.

cell_features = cell(img_n,1);  cell_validpts = cell(img_n,1); 
orig_lines = cell(img_n,1);
fprintf('  Keypoint detection and matching...');tic;
for i=1:img_n
    tmppoints = detectSIFTFeatures(rgb2gray(cell_imgs{i}));
    [cell_features{i}, cell_validpts{i}] = extractFeatures(rgb2gray(cell_imgs{i}),tmppoints);
    orig_lines{i}=lsd(cell_paths{i});
end


% find matches between different image pair
inliers = cell(size(imgs_pairs,1),1);
data_norm = cell(size(imgs_pairs,1),1);
data_orig = cell(size(imgs_pairs,1),1);
data_lines = cell(size(imgs_pairs,1),1);
h  = cell(size(imgs_pairs,1),1);
for k=1:size(imgs_pairs,1)
    m = imgs_pairs(k,1);  n = imgs_pairs(k,2);
    validpts1=cell_validpts{m}; validpts2=cell_validpts{n};
    indexPairs = matchFeatures(cell_features{m}, cell_features{n}, 'MatchThreshold',10);
    matchedPoints1 = validpts1(indexPairs(:,1),:);
    matchedPoints2 = validpts2(indexPairs(:,2),:);
    pts1=double(matchedPoints1.Location');
    pts2=double(matchedPoints2.Location');
    data_orig{k} = [ pts1 ; ones(1,size(pts1,2)) ; pts2 ; ones(1,size(pts2,2)) ]; 
    
    %% APAP's multi-sampling method
    % outlier removal - Multi-GS (RANSAC).
    % normalise point distribution.
    [data_norm_img1, ~] = normalise2dpts(data_orig{k}(1:3,:));
    [data_norm_img2, ~] = normalise2dpts(data_orig{k}(4:6,:));
    data_norm{k} = [data_norm_img1; data_norm_img2];
    rng(0);
    [ ~,res,~,~, err ] = multigsSampling(100, data_norm{k}, 500, 10);
    con = sum(res<=0.05);
    [ ~, maxinx ] = max(con);
    inliers{k} = find(res(:,maxinx)<=0.05);

    % match lines based on point inliers
    [tmpline1, tmpline2] = multiLineMatch(cell_imgs{m},cell_imgs{n},orig_lines{m},orig_lines{n},...
        data_orig{k}(1:2,inliers{k}),data_orig{k}(4:5,inliers{k})); 
   
    data_lines{k} = [tmpline1, tmpline2]; 
    
    if isempty(tmpline1)
        h{k} = calcHomo(data_orig{k}(1:2,inliers{k}),data_orig{k}(4:5,inliers{k}));
        h{k} = h{k}(:);
    else
        % refine homography using point inliers and line matches.
        h{k} = calcHomoPointLine(data_orig{k}(1:2,inliers{k}),data_orig{k}(4:5,inliers{k}),...
                                tmpline1, tmpline2);
        h{k} = h{k}(:); % <- these are the initial estimates for Levenberg-Marquardt.
    end
%     % verify image match.
%     if err || length(inliers{k}) <= 5.9 + 0.22 * length(matches)%(8 + 0.3 *length(matches)) % invalid image match (following Brown and Lowe).
%                                                       % conference paper version values: 5.9 + 0.22 * length(matches)
%         % mark image pair as invalid image match and clean up the image pair's data.
%         imgs_pairs(k,3) = 0;
%     end
    
end
fprintf('done (%fs)\n',toc);

end

