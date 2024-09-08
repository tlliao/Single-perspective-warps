clear; clc; close all;
%% Setup VLFeat toolbox and add other papers' codes.
%----------------------
addpath('modelspecific'); addpath('multigs');
addpath('texture_mapping');
addpath('LSD_matlab'); addpath('LineMatching');
run('vlfeat-0.9.21/toolbox/vl_setup');
grid_size = 40;
%----------------- parameter of energy minimization
%     line-align-(al), shape-(ps), shorten-(pj), line-(le)
parameters = [5, 50, 5, 5];
%------------------------------
%------------------
% Multiple images to stitch.
%------------------
%% Initialization
pathname = strcat('Imgs\...\');
outpath = strcat(pathname, 'results\'); %-testparas
imgs_format = '*.jpg';
dir_folder = dir(strcat(pathname, imgs_format));
img_n = length(dir_folder);  % number of files (with specified format) in the folder.
if img_n<=2; error('Input images are less than 3! Please add more images.'); end 
if ~exist(outpath,'dir'); mkdir(outpath); end   
%-------------

%% Read images.
%-------------
fprintf('> Reading images...');tic;
cell_imgs = cell(img_n, 1);  cell_paths = cell(img_n,1);
for i=1:img_n
    cell_paths{i} = sprintf('%s%s',pathname,dir_folder(i).name);
    cell_imgs{i}  = im2double(imread(cell_paths{i}));
end
fprintf('done (%fs)\n',toc);

% Resolution/grid-size for the mapping function (divide it into C1*C2 cells).
C1 = ceil(size(cell_imgs{img_n},1)/grid_size);
C2 = ceil(size(cell_imgs{img_n},2)/grid_size);
num_V = (C1+1)*(C2+1); % number of control vertices in each image

%% detect and match sift features
[data_orig, inliers, data_lines, h, imgs_pairs] = pointlineMatch(cell_imgs, cell_paths); 

%% our single-perspective warp and blending 
fprintf('  our SPW warp and blending...');tic;
[ baH, ref_n, points_xk, points_predict, lines_xk, lines_predict] = bundleAdjustmentRefofLine(data_orig, inliers, data_lines, h, imgs_pairs);

%% generate u and v sample points of all target images (img_n-1 targets)
[ cell_us, cell_ue, cell_v, nor_vecs ] = generateUVsamples( cell_imgs, C1, C2, baH, ref_n );

%% calculate the v energy term (v sample points: slope-preserving and part equidistant)
%fprintf('> generating v sample points and v term to optimize V*...');tic;
sparse_v = energyLineV(cell_imgs, C1, C2, cell_v, nor_vecs);
%fprintf('done (%fs)\n',toc);

%% calculate the u energy term (u sample points: slope-preserving and part equidistant)
%fprintf('> generating u sample points and u term to optimize V*...');tic;
[sparse_us, sparse_ue] = energyLineU(cell_imgs, C1, C2, baH, cell_us, cell_ue );
%fprintf('done (%fs)\n',toc);

%% calculate the new alignment energy term with scale operator
%fprintf('> generating scale-alignment term ||AV*-b||^2 to optimize V*...');tic;
[ sparse_point_align, sparse_line_align, pMatch_ones, cMatch_ones, initial_pc ] = energyPointLineAlign(cell_imgs, points_xk, points_predict, lines_xk, lines_predict, C1, C2, ref_n);


%% calculate line-preserving term with line segments
%fprintf('> detect line segments and calculate line-preserving term...');tic;
[cell_lines, cell_slope] = linesDetect( cell_paths, cell_imgs, C1, C2, ref_n );  % detect lines using LSD method
sparse_line = energyLineSegment(cell_imgs, cell_lines, cell_slope, baH, C1, C2);
%fprintf('done (%fs)\n', toc);

%% construct matrix A,b to minimize ||Ax-b||^2  (A'*Ax=A'*b)
zero_len = size(sparse_us,1) + size(sparse_v,1) + size(sparse_ue,1)+ size(sparse_line,1);

%     V0 = zeros(2*num_V*(img_n-1),1); 
% initilize gauss: Homography
V0=zeros(2*num_V*(img_n-1),1);k=1;
for i=1:img_n
    if i==ref_n; continue; end
    [X, Y] = meshgrid(linspace(1,size(cell_imgs{i},2),C2+1), linspace(1,size(cell_imgs{i},1),C1+1)); % mesh grid index
    Mv = [X(:), Y(:)];
    tmpH = baH{ref_n}/baH{i};
    warp_Mv = tmpH*[Mv'; ones(1,length(Mv))]; warp_Mv=warp_Mv(1:2,:)./repmat(warp_Mv(3,:),2,1);
    V0((k-1)*2*num_V+1:k*2*num_V)=warp_Mv(:);
    k=k+1;
end

len_ones = size(sparse_point_align,2)-2*(img_n-1)*num_V;  % number of p and c to be optimized

for parai=1:size(parameters,1)
    para_l=parameters(parai,1); para_ps=parameters(parai,2); 
    para_pj=parameters(parai,3); para_s=parameters(parai,4);
    solve_A = [sparse_point_align;
               sqrt(para_l).*sparse_line_align;
               sqrt(para_ps).*sparse_us, zeros(size(sparse_us,1),len_ones); 
               sqrt(para_ps).*sparse_v, zeros(size(sparse_v,1), len_ones); 
               sqrt(para_pj).*sparse_ue, zeros(size(sparse_ue,1),len_ones);
               sqrt(para_s).*sparse_line, zeros(size(sparse_line,1),len_ones)];
    m_x = [pMatch_ones; sqrt(para_l).*cMatch_ones; zeros(zero_len,1)]; 

   %% use iterative methods in sparse linear system to solve the energy minimization        
    fprintf('  Use LSQR method to calcuate optimized V_star...');tic;
    [V_star, flag, ~, iter, resvec] = lsqr(solve_A, m_x, 1e-8, 10000, [],[], [V0; initial_pc]);
    fprintf('done (%fs)\n',toc);
    optimized_V = vec2mat(V_star(1:2*(img_n-1)*num_V), 2);  clear V_star;

   %% construt a canvas to map image on
    % % Output Canvas
    off = round([ 1 - min([1 optimized_V(:,1)']) + 1 ; 1 - min([1 optimized_V(:,2)']) + 1 ]);
    cw = max( size(cell_imgs{ref_n},2), ceil(max(optimized_V(:,1))) ) + off(1)-1 ;
    ch = max( size(cell_imgs{ref_n},1), ceil(max(optimized_V(:,2))) ) + off(2)-1 ;
    if cw>1e4 || ch>1e4
        continue;%output_canvas = cell_imgs{ref};  % if the warp is too sharp, return the reference image immediately
    else
        cell_wimgs = cell(img_n-1,1); cell_Masks = cell(img_n-1,1);
        for i=1:img_n-1
            tmp_opV = optimized_V(num_V*(i-1)+1:num_V*i, :); % the warped control vertices in axis
            wX1 = reshape(tmp_opV(:,1), C1+1, C2+1);  wY1 = reshape(tmp_opV(:,2), C1+1, C2+1);
            if i<ref_n
                ori_img = cell_imgs{i};
            else 
                ori_img = cell_imgs{i+1};
            end
            [X1, Y1] = meshgrid(linspace(1,size(ori_img,2),C2+1), linspace(1,size(ori_img,1),C1+1));
            cell_wimgs{i} = meshmap_warp2homo(ori_img, X1, Y1, wX1, wY1);
            cell_Masks{i} = imfill(imbinarize(rgb2gray(cell_wimgs{i}), 0),'holes');
        end
        cell_Homo = cell(img_n,1); cell_HomoMask = cell(img_n,1);
        k=1;
        for i=1:img_n
           cell_Homo{i} = zeros(ch,cw,3); cell_HomoMask{i} = false(ch,cw);
           if i==ref_n
               cell_Homo{i}(off(2):off(2)-1+size(cell_imgs{ref_n},1), off(1):off(1)-1+size(cell_imgs{ref_n},2),:) = cell_imgs{ref_n};
               cell_HomoMask{i}(off(2):off(2)-1+size(cell_imgs{ref_n},1), off(1):off(1)-1+size(cell_imgs{ref_n},2)) = true;
           else
               tmp_V = optimized_V(num_V*(k-1)+1:num_V*k, :);
               cell_Homo{i}(round(min(tmp_V(:,2)))+off(2)-1:round(min(tmp_V(:,2)))+off(2)-2+size(cell_wimgs{k},1),...
                            round(min(tmp_V(:,1)))+off(1)-1:round(min(tmp_V(:,1)))+off(1)-2+size(cell_wimgs{k},2), :) = cell_wimgs{k};
               cell_HomoMask{i}(round(min(tmp_V(:,2)))+off(2)-1:round(min(tmp_V(:,2)))+off(2)-2+size(cell_wimgs{k},1),...
                            round(min(tmp_V(:,1)))+off(1)-1:round(min(tmp_V(:,1)))+off(1)-2+size(cell_wimgs{k},2)) = cell_Masks{k};
                k=k+1;
           end
        end
        output = imageBlending(cell_Homo, cell_HomoMask, 'linear');
        clear cell_Homo   cell_HomoMask;
    end
    pngpara=sprintf('%d-%d-%d-%d',parameters(parai,1),parameters(parai,2),parameters(parai,3),parameters(parai,4));
    imwrite(output, [outpath,'linear-',pngpara,'.png']); 
end
fprintf('done (%fs)',toc);