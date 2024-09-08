function [ cell_us, cell_ue, cell_v, nor_vecs ] = generateUVsamples( cell_imgs, C1, C2, baH, ref_n )
% for each target image, generate u,v sample points based on global homography (baH) r.s.t 
% reference image after bundle adjustment
% the u's equidistant sampling is based on global homography baH
% this version is designed for a sequence of input images
num_imgs = size(baH,1);  % number of total images (include reference)
cell_us = cell(num_imgs, 1);    % u-slope sample points  
cell_ue = cell(num_imgs, 1);    % u-equal sample points
cell_v =  cell(num_imgs, 1);    % v-slope-equal sample points
nor_vecs = zeros(num_imgs,2);   % normal vectors of v sample points
multi = 2;                      % density of sample points r.t. mesh grids 

%% generating u,v sample points for all targets
for numi=[1:ref_n-1, ref_n+1:num_imgs]
    init_H = baH{ref_n}/baH{numi};  init_H = init_H./init_H(end);
    nor_vecs(numi,:) = [init_H(2)*init_H(6)-init_H(5)*init_H(3), init_H(4)*init_H(3)-init_H(1)*init_H(6)];
    nor_vecs(numi,:) = nor_vecs(numi,:)./norm(nor_vecs(numi,:));
    theta = atan2(-init_H(6),-init_H(3));
    [M, N, ~] = size(cell_imgs{numi});
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    off_center = [(N+1)/2; (M+1)/2] - R*[(N+1)/2; (M+1)/2];
    [X, Y] = meshgrid(linspace(2-N, 2*N-1, 3*multi*(C2+1)-2), linspace(2-M, 2*M-1, 3*multi*(C1+1)-2)); % mesh grid index
    
    lines_v = zeros(2*size(X,2), size(X,1));
    k=1;
    % rotated vertical line slope-preserving
    for j=1:size(X,2)
        tmp_Mv = [X(:,j)'; Y(:,j)'];  % each vertical line
        tmp_line = R*tmp_Mv + repmat(off_center, 1,length(tmp_Mv)); % rotated vertical line
        inner_ind = (tmp_line(1,:)>=1) & (tmp_line(1,:)<=N) & (tmp_line(2,:)>=1) & (tmp_line(2,:)<=M); % rotated vertical line lies in img
        lines_v(2*k-1:2*k, 1:sum(inner_ind)) = tmp_line(:,inner_ind);  % useful rotated vertical line
        lines_v(2*k-1:2*k, end) = sum(inner_ind);  % number of sample points on the line
        k=k+1;
    end 
    lines_v(all(sum(lines_v,2),2)==0,:)=[];  % delete all-zero rows
    lines_v(:,all(sum(lines_v,1),1)==0)=[];  % delete all-zero columns
    cell_v{numi} = lines_v; % left 1-neighbor target v sample points

    lines_us = zeros(2*size(X,1), size(X,2));
    k=1;
    % rotated vertical line slope-preserving
    for i=1:size(X,1)
        tmp_Mv = [X(i,:); Y(i,:)];  % each vertical line
        tmp_line = R*tmp_Mv + repmat(off_center, 1, length(tmp_Mv)); % rotated vertical line
        inner_ind = (tmp_line(1,:)>=1) & (tmp_line(1,:)<=N) & (tmp_line(2,:)>=1) & (tmp_line(2,:)<=M); % rotated vertical line lies in img
        lines_us(2*k-1:2*k, 1:sum(inner_ind)) = tmp_line(:,inner_ind);  % useful rotated vertical line
        lines_us(2*k-1:2*k, end) = sum(inner_ind);  % number of sample points on the line
        k=k+1;
    end
    lines_us(all(sum(lines_us,2),2)==0,:)=[];  % delete all-zero rows
    lines_us(:,all(sum(lines_us,1),1)==0)=[];  % delete all-zero columns
    cell_us{numi} = lines_us; % left 1-neighbor target u-slope sample points

   %% filter u sample points (omit points in the overlapping region) based on baH
    if numi<=ref_n-1
        neigh_H = baH{ref_n}/baH{numi+1};
        sz1 = size(cell_imgs{numi+1},1);
        neigh_pts1 = neigh_H*[1;-sz1;1];   neigh_pts1 = neigh_pts1(1:2)./neigh_pts1(3);
        neigh_pts2 = neigh_H*[1;2*sz1;1];  neigh_pts2 = neigh_pts2(1:2)./neigh_pts2(3);
        left_or_right = 1;
    else
        neigh_H = baH{ref_n}/baH{numi-1};
        [sz1,sz2,~] = size(cell_imgs{numi-1});
        neigh_pts1 = neigh_H*[sz2;-sz1;1];   neigh_pts1 = neigh_pts1(1:2)./neigh_pts1(3);
        neigh_pts2 = neigh_H*[sz2;2*sz1;1];  neigh_pts2 = neigh_pts2(1:2)./neigh_pts2(3);
        left_or_right = 0;        
    end

    x1=neigh_pts1(1); y1=neigh_pts1(2);  x2=neigh_pts2(1); y2=neigh_pts2(2);
    lines_ue = zeros(size(lines_us));  % u sample points to preserve equidistant
    for i=1:2:size(lines_us,1)-1
        num_v = lines_us(i,end);
        warp_us = init_H*[lines_us(i:i+1,1:num_v);ones(1,num_v)]; warp_us = warp_us(1:2,:)./repmat(warp_us(3,:),2,1);
        vec_prod = (x1-warp_us(1,:)).*(y2-warp_us(2,:))-(y1-warp_us(2,:)).*(x2-warp_us(1,:));
        left = vec_prod>0;  right = vec_prod<0;
        index_prod = (left_or_right & left) | ((~left_or_right) & right);  % find the sample points in the non-overlapping region
        lines_ue(i:i+1, 1:sum(index_prod)) = lines_us(i:i+1, index_prod);
        lines_ue(i:i+1, end) = sum(index_prod); % number of sample points
    end
    lines_ue(all(sum(lines_ue,2),2)==0,:)=[];  % delete all-zero rows
    lines_ue(:,all(sum(lines_ue,1),1)==0)=[];  % delete all-zero columns
    if isempty(lines_ue); lines_ue=[0;0]; end
    cell_ue{numi} = lines_ue;   % left 1-neighbor target u-equal sample points
end


end