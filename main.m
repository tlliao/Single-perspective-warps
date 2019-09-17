clear; clc; close all;
%----------------------
addNeedPaths; % add paths and toolbox
%--------------------------

% Parameters of energy minimization (mesh deformation)
parameters.grid_size = 40;  % grid size of mesh deformation

parameters.line_align = 5;  % line alignment term
parameters.perspective = 50; % perspective term
parameters.projective = 5;   % projective term
parameters.saliency = 5;    % saliency term

parameters.line_threshold = 50;  % threshold for length of matches lines
%------------------------
% Images to stitch.
%-----------------------
tar = 2; ref_n = 1;  % choose the target and reference image
pathname = strcat('Imgs\'); %
outpath = strcat(pathname, 'results\');
imgs_format = '*.jpg'; 
dir_folder = dir(strcat(pathname, imgs_format));
if ~exist(outpath,'dir'); mkdir(outpath); end
path1 =  sprintf('%s%s',pathname,dir_folder(tar).name); %
path2 =  sprintf('%s%s',pathname,dir_folder(ref_n).name); %
%-------------
%% Read images.
%-------------
fprintf('> Reading images...');tic;
img1 = im2double(imread(path1));   img2 = im2double(imread(path2));
fprintf('done (%fs)\n',toc);
% Resolution/grid-size for the mapping function (divide it into C1*C2 cells).
C1 = ceil(size(img1,1)/parameters.grid_size);
C2 = ceil(size(img1,2)/parameters.grid_size);

%% detect and match sift features
[pts1, pts2] = siftMatch(img1, img2); 
[matches_1, matches_2] = multiSample_APAP(pts1, pts2);

%% detect and match line segments
[line_match1, line_match2] = twoLineMatch(path1, path2, matches_1, matches_2, parameters);

%% our single-perspective warp (SPW) and blending
fprintf('  Our SPW warp and blending...');tic;
[h, ~, T1, T2] = calcHomoPointLine( matches_1, matches_2, line_match1, line_match2 );
pts_line_H = T2\(h*T1);

%% generating mesh grid (C1*C2) to optimize warped control vertices and rotation angle theta
[X, Y] = meshgrid(linspace(1,size(img1,2),C2+1), linspace(1,size(img1,1),C1+1)); % mesh grid index
% Mesh (cells) vertices' coordinates.
Mv = [X(:), Y(:)];
init_H = pts_line_H;
init_H = init_H./(init_H(end));
theta = atan2(-init_H(6), -init_H(3));  % evaluate the rotation (r.t. original coordinate system)

%% use rotation angle theta to calculate normal vector of warped v-line normal_vec
fprintf('> generating u-v sample points and u-v term to optimize V*...');tic;
[lines_vs, lines_us, lines_ue] = generateUV( img1, img2, init_H, theta, C1 ,C2 ); % rotated vertical and horizontal lines
nor_vec_v = [init_H(2)*init_H(6)-init_H(5)*init_H(3), init_H(4)*init_H(3)-init_H(1)*init_H(6)]; % the normal vector of v-lines after transformation 
nor_vec_v = nor_vec_v./norm(nor_vec_v);  % normalization
sparse_v  = energyLineV( img1, C1, C2, lines_vs, nor_vec_v ); % energy of preserving v-lines
[sparse_us, sparse_ue] = energyLineU( img1, C1, C2, lines_us, lines_ue, init_H );  % energy of preserving u-slope and u-equidistant
fprintf('done (%fs)\n',toc);

%% calculate the new alignment energy term with scale operator
fprintf('> generating scale-alignment term ||AV*-b||^2 to optimize V*...');tic;
[ sparse_align, psMatch ] = energyAlign( img1, C1, C2, matches_1, matches_2 );
[ sparse_line_align, cMatch ] = energyLineAlign( img1, C1, C2, line_match1, line_match2 );
fprintf('done (%fs)\n',toc);

%% calculate line-preserving term with line segments
fprintf('> detect line segments and calculate line-preserving term to optimize V*...');tic;
[sa_lines, sl_lines] = linesDetect( path1, img1, C1, C2 );  % detect lines using LSD method
sparse_line = energyLineSegment(img1, sa_lines, sl_lines, init_H, C1, C2);
fprintf('done (%fs)\n', toc);

%% construct matrix A,b to minimize ||Ax-b||^2  (A'*Ax=A'*b)
%  E_warp = E_align + E_usample + E_vsample + E_line, 
zero_len = size(sparse_us,1)+size(sparse_v,1)+size(sparse_ue,1)+size(sparse_line,1);   

% generating initial warping hypothesis (homography)
warp_hv = init_H*[Mv'; ones(1,length(Mv))];
warp_hv = warp_hv(1:2,:)./repmat(warp_hv(3,:),2,1);
init_V = warp_hv(:);

para_l = parameters.line_align; para_ps = parameters.perspective;
para_pj = parameters.projective; para_s = parameters.saliency; 
Matrix_A = [sparse_align; sqrt(para_l).*sparse_line_align; sqrt(para_ps).*sparse_us; sqrt(para_ps).*sparse_v;...
            sqrt(para_pj).*sparse_ue;  sqrt(para_s).*sparse_line]; 
m_x = [psMatch; sqrt(para_l).*cMatch; zeros(zero_len,1)]; 
%% use iterative methods in sparse linear system to solve the energy minimization
fprintf('> Use LSQR method to calcuate optimized V_star...');tic;
[V_star, flag, ~, iter] = lsqr(Matrix_A, m_x, 1e-8 , 5000, [], [], init_V);
fprintf('done (%fs)\n',toc);
optimized_V = vec2mat(V_star,2); clear V_star; % show the warped control vertices in axis

%% calculate the warp image using warp function (homography)
fprintf('> mesh deformation using bilinear interpolation...');tic;
wX = reshape(optimized_V(:,1), C1+1, C2+1);
wY = reshape(optimized_V(:,2), C1+1, C2+1);
warped_img1 = meshmap_warp2homo(img1, X, Y, wX, wY);
fprintf('done (%fs)\n',toc);

%%  blend the warped images, show the final result (compared to Homography)
% % Canvas size.
off = ceil([ 1 - min([1 optimized_V(:,1)']) + 1 ; 1 - min([1 optimized_V(:,2)']) + 1 ]);
cw = max([ceil(optimized_V(:,1))', size(img2,2)])+off(1)-1;
ch = max([ceil(optimized_V(:,2))', size(img2,1)])+off(2)-1;

img1Homo = zeros(ch,cw,3); img2Homo = zeros(ch,cw,3);
img1Homo(floor(min(optimized_V(:,2)))+off(2)-1:floor(min(optimized_V(:,2)))+off(2)-2+size(warped_img1,1),...
    floor(min(optimized_V(:,1)))+off(1)-1:floor(min(optimized_V(:,1)))+off(1)-2+size(warped_img1,2), :) = warped_img1; 
img2Homo(off(2):(off(2)+size(img2,1)-1),off(1):(off(1)+size(img2,2)-1),:) = img2;

linear_out = imageBlending(img1Homo, img2Homo, 'linear');       
pngout = sprintf('linear-%d-%d-%d-%d.jpg', para_l, para_ps, para_pj, para_s); 
imwrite(linear_out, [outpath, pngout]);
    
fprintf('done (%fs)\n',toc);

