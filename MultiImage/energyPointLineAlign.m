function [ sparse_point_align, sparse_line_align, pMatch_ones, cMatch_ones, initial_pc ] = energyPointLineAlign(cell_imgs, points_xk, points_predict, lines_xk, lines_predict, C1, C2, ref_n)
%   for given target image, generate C1*C2 mesh grid for points align 
%   sparse_point_align is a sparse matrix A to solve Ax-b=0
%   points_xk: 2*num_imgs x M matrix
%   points_predict: 2 x M matrix 
%   M is the number of point predictions
%   we do not need to delete the all-zero rows  in A and crossponding b, as
%   A'*A is equal w or w/o all-zero

%% average alignment: with 1/sum(\delta_ik)
num_imgs = size(cell_imgs,1);  num_V = (C1+1)*(C2+1);
num_points = size(points_xk, 2);

indexofpoints = ~isnan(points_xk);             
indexofpoints = indexofpoints(2:2:end,:);   % index of all images in points_xk

X_col = linspace(1, size(cell_imgs{ref_n},2), C2+1); % column index of cells
Y_row = linspace(1, size(cell_imgs{ref_n},1), C1+1); % row index of cells
x_dis = X_col(2)-X_col(1);  % the width of cell
y_dis = Y_row(2)-Y_row(1);  % the height of cell

% align all target points to reference frame (including matches w/ and w/o reference image)
k=1; kk=1; js=1;
Mesh_ps = zeros(4,2);
sp_i=zeros(8*num_imgs*num_points, 1);
sp_j=sp_i; sp_s = sp_i;
pMatch = zeros(2*num_imgs*num_points, 1);
sp_ii=zeros(8*num_imgs*num_points, 1);
sp_jj=sp_ii; sp_ss = sp_ii;
sp_pi = zeros(4*num_imgs*num_points, 1); sp_pj = sp_pi; sp_ps = sp_pi;
initial_p = zeros(2*num_points, 1);  % nitial gauss of c
for i=1:num_points
    if indexofpoints(ref_n, i)  % matches w/ reference
       tmpindex = setdiff(find(indexofpoints(:,i)), ref_n); 
       ref_px = points_predict(1,i); ref_py = points_predict(2,i);
       for j=tmpindex'
           index_j=j; if j>ref_n; index_j=j-1; end
           px1 = points_xk(2*j-1, i); py1 = points_xk(2*j, i);
           cx1 = min(find( (px1-X_col)<x_dis  & (px1-X_col)>=0, 1), C2); % the x index of endpoint's position
           cy1 = min(find( (py1-Y_row)<y_dis  & (py1-Y_row)>=0, 1), C1); % the y index of endpoint's position
           % the cell containing point p
           Mesh_ps(1:4,:) = [X_col(cx1), Y_row(cy1);     % v1
                               X_col(cx1+1), Y_row(cy1);   % v2
                               X_col(cx1+1), Y_row(cy1+1); % v3
                               X_col(cx1), Y_row(cy1+1)];   % v4
           coeff_mesh_ps = meshGridAlign(Mesh_ps, [px1,py1]);
           nums1 = (C1+1)*(cx1-1)+cy1; % index of v1*       
           nums2 = nums1+(C1+1);     % index of v2*        
           nums3 = nums2+1;            % index of v3*
           nums4 = nums1+1;            % index of v4*
           sp_i(8*k-7: 8*k) = [ones(1,4).*(2*k-1), ones(1,4).*(2*k)];
           sp_j(8*k-7: 8*k) = [2*nums1-1, 2*nums2-1, 2*nums3-1, 2*nums4-1,...
               2*nums1, 2*nums2, 2*nums3, 2*nums4] + 2*num_V*(index_j-1);
           sp_s(8*k-7: 8*k) = 1/sqrt(length(tmpindex)).*[coeff_mesh_ps; coeff_mesh_ps];
           pMatch(2*k-1:2*k) = 1/sqrt(length(tmpindex)).*[ref_px; ref_py]; 
           k=k+1;
       end
    else   % matches w/o reference
       tmpindex = find(indexofpoints(:,i)); 
       pre_px = points_predict(1,i); pre_py = points_predict(2,i);
       for j=tmpindex'
           index_j=j; if j>ref_n; index_j=j-1; end
           px1=points_xk(2*j-1, i); py1=points_xk(2*j, i);
           cx1 = min(find( (px1-X_col)<x_dis  & (px1-X_col)>=0, 1), C2); % the x index of endpoint's position
           cy1 = min(find( (py1-Y_row)<y_dis  & (py1-Y_row)>=0, 1), C1); % the y index of endpoint's position
           % the cell containing point p
           Mesh_ps(1:4,:) = [X_col(cx1), Y_row(cy1);     % v1
                               X_col(cx1+1), Y_row(cy1);   % v2
                               X_col(cx1+1), Y_row(cy1+1); % v3
                               X_col(cx1), Y_row(cy1+1)];   % v4
           coeff_mesh_ps = meshGridAlign(Mesh_ps, [px1,py1]);
           nums1 = (C1+1)*(cx1-1)+cy1; % index of v1*       
           nums2 = nums1+(C1+1);     % index of v2*        
           nums3 = nums2+1;            % index of v3*
           nums4 = nums1+1;            % index of v4*
           sp_ii(8*kk-7: 8*kk) = [ones(1,4).*(2*kk-1), ones(1,4).*(2*kk)];
           sp_jj(8*kk-7: 8*kk) = [2*nums1-1, 2*nums2-1, 2*nums3-1, 2*nums4-1,...
               2*nums1, 2*nums2, 2*nums3, 2*nums4] + 2*num_V*(index_j-1);
           sp_ss(8*kk-7: 8*kk) = 1/sqrt(length(tmpindex)).*[coeff_mesh_ps; coeff_mesh_ps];
           
           sp_pi(4*kk-3:4*kk) = [2*kk-1,2*kk-1,2*kk,2*kk];
           sp_pj(4*kk-3:4*kk) = [2*js-1,2*js,2*js-1,2*js];
           sp_ps(4*kk-3:4*kk) = [-1,0;0,-1];
           
           kk=kk+1;     
       end
       initial_p(2*js-1:2*js) = 1/sqrt(length(tmpindex)).*[pre_px; pre_py];
       js=js+1;
    end
end

% delete redundant zeros in sp_i/sp_ii, sp_j/sp_jj, sp_s/sp_ss, c_Match and c_ones
sp_i = sp_i(1:8*(k-1)); sp_j = sp_j(1:8*(k-1)); sp_s = sp_s(1:8*(k-1));
pMatch = pMatch(1:2*(k-1));
sp_ii = sp_ii(1:8*(kk-1)); sp_jj = sp_jj(1:8*(kk-1)); sp_ss = sp_ss(1:8*(kk-1));
sp_pi = sp_pi(1:4*(kk-1)); sp_pj = sp_pj(1:4*(kk-1)); sp_ps= sp_ps(1:4*(kk-1));
initial_p = initial_p(1:2*(js-1));
sparse_point_ref = sparse(sp_i, sp_j, sp_s, max(sp_i), 2*(num_imgs-1)*num_V);
if isempty(sp_ii)
    p_ones = [];
    sparse_point_free=[];
else
    p_ones = sparse(sp_pi, sp_pj, sp_ps, max(sp_pi), 2*(js-1)); 
    sparse_point_free = sparse(sp_ii, sp_jj, sp_ss, max(sp_ii), 2*(num_imgs-1)*num_V);
end

% % final sparse matrix of point alignment
% sparse_point_align = [sparse_point_ref, zeros(max(sp_i), size(p_ones,2));
%     sparse_point_free, p_ones];
% pMatch_ones = [pMatch; zeros(size(p_ones, 1),1)];

%% average alignment: with 1/sum(\delta_ik)
num_lines = size(lines_xk, 2);

indexoflines = ~isnan(lines_xk);             
indexoflines = indexoflines(4:4:end,:);   % index of all images in lines_xk

abc_predict = [lines_predict(4,:)-lines_predict(2,:); lines_predict(1,:)-lines_predict(3,:);
    lines_predict(3,:).*lines_predict(2,:)-lines_predict(1,:).*lines_predict(4,:)];
abc_predict = abc_predict./repmat(sqrt(abc_predict(1,:).^2+abc_predict(2,:).^2), 3, 1);

% align all target lines to reference frame (including matches w/ and w/o reference image)
k=1; kk=1; js=1;
Mesh_ps = zeros(4,2); Mesh_pe = zeros(4,2);
sp_i=zeros(16*num_imgs*num_lines, 1);
sp_j=sp_i; sp_s = zeros(16*num_imgs*num_lines, 1);
cMatch = zeros(2*num_imgs*num_lines, 1);
sp_ii=zeros(16*num_imgs*num_lines, 1);
sp_jj=sp_ii; sp_ss = zeros(16*num_imgs*num_lines, 1);
sp_ci = zeros(2*num_imgs*num_lines, 1); sp_cj=sp_ci; sp_cs = sp_ci;
initial_c = zeros(num_lines, 1);  % nitial gauss of c
for i=1:num_lines
    if indexoflines(ref_n,i)  % matches w/ reference
       tmpindex = setdiff(find(indexoflines(:,i)),ref_n); 
       a=abc_predict(1,i); b=abc_predict(2,i); c=abc_predict(3,i);
       for j=tmpindex'
           index_j=j; if j>ref_n; index_j=j-1; end
           px1=lines_xk(4*j-3,i); py1=lines_xk(4*j-2,i);
           px2=lines_xk(4*j-1,i); py2=lines_xk(4*j,i);
           cx1 = min(find( (px1-X_col)<x_dis  & (px1-X_col)>=0, 1), C2); % the x index of endpoint's position
           cy1 = min(find( (py1-Y_row)<y_dis  & (py1-Y_row)>=0, 1), C1); % the y index of endpoint's position
           cx2 = min(find( (px2-X_col)<x_dis  & (px2-X_col)>=0, 1), C2); % the x index of endpoint's position
           cy2 = min(find( (py2-Y_row)<y_dis  & (py2-Y_row)>=0, 1), C1); % the y index of endpoint's position
           % the cell containing line segments p
           Mesh_ps(1:4,:) = [X_col(cx1), Y_row(cy1);     % v1
                               X_col(cx1+1), Y_row(cy1);   % v2
                               X_col(cx1+1), Y_row(cy1+1); % v3
                               X_col(cx1), Y_row(cy1+1)];   % v4
           Mesh_pe(1:4,:) = [X_col(cx2), Y_row(cy2);     % v1
                               X_col(cx2+1), Y_row(cy2);   % v2
                               X_col(cx2+1), Y_row(cy2+1); % v3
                               X_col(cx2), Y_row(cy2+1)];   % v4
           coeff_mesh_ps = meshGridAlign(Mesh_ps, [px1,py1]);
           coeff_mesh_pe = meshGridAlign(Mesh_pe, [px2,py2]);
           nums1 = (C1+1)*(cx1-1)+cy1; % index of v1*       
           nums2 = nums1+(C1+1);     % index of v2*        
           nums3 = nums2+1;            % index of v3*
           nums4 = nums1+1;            % index of v4*
           nume1 = (C1+1)*(cx2-1)+cy2; % index of v1*       
           nume2 = nume1+(C1+1);     % index of v2*        
           nume3 = nume2+1;            % index of v3*
           nume4 = nume1+1;            % index of v4*  
           sp_i(16*k-15:16*k)=[ones(1,8).*(2*k-1), ones(1,8).*(2*k)];
           sp_j(16*k-15:16*k)=[2*nums1-1, 2*nums2-1, 2*nums3-1, 2*nums4-1,...
               2*nums1, 2*nums2, 2*nums3, 2*nums4,...
               2*nume1-1, 2*nume2-1, 2*nume3-1, 2*nume4-1,...
               2*nume1, 2*nume2, 2*nume3, 2*nume4] + 2*num_V*(index_j-1);
           sp_s(16*k-15:16*k)= 1/sqrt(length(tmpindex)).*[a.*coeff_mesh_ps; b.*coeff_mesh_ps; a.*coeff_mesh_pe; b.*coeff_mesh_pe];
           cMatch(2*k-1:2*k) = 1/sqrt(length(tmpindex)).*[-c;-c]; 
           k=k+1;
       end
    else   % matches w/o reference
       tmpindex = find(indexoflines(:,i)); 
       a=abc_predict(1,i); b=abc_predict(2,i); c=abc_predict(3,i);
       for j=tmpindex'
           index_j=j; if j>ref_n; index_j=j-1; end
           px1=lines_xk(4*j-3,i); py1=lines_xk(4*j-2,i);
           px2=lines_xk(4*j-1,i); py2=lines_xk(4*j,i);
           cx1 = min(find( (px1-X_col)<x_dis  & (px1-X_col)>=0, 1), C2); % the x index of endpoint's position
           cy1 = min(find( (py1-Y_row)<y_dis  & (py1-Y_row)>=0, 1), C1); % the y index of endpoint's position
           cx2 = min(find( (px2-X_col)<x_dis  & (px2-X_col)>=0, 1), C2); % the x index of endpoint's position
           cy2 = min(find( (py2-Y_row)<y_dis  & (py2-Y_row)>=0, 1), C1); % the y index of endpoint's position
           % the cell containing line segments p
           Mesh_ps(1:4,:) = [X_col(cx1), Y_row(cy1);     % v1
                               X_col(cx1+1), Y_row(cy1);   % v2
                               X_col(cx1+1), Y_row(cy1+1); % v3
                               X_col(cx1), Y_row(cy1+1)];   % v4
           Mesh_pe(1:4,:) = [X_col(cx2), Y_row(cy2);     % v1
                               X_col(cx2+1), Y_row(cy2);   % v2
                               X_col(cx2+1), Y_row(cy2+1); % v3
                               X_col(cx2), Y_row(cy2+1)];   % v4
           coeff_mesh_ps = meshGridAlign(Mesh_ps, [px1,py1]);
           coeff_mesh_pe = meshGridAlign(Mesh_pe, [px2,py2]);
           nums1 = (C1+1)*(cx1-1)+cy1; % index of v1*       
           nums2 = nums1+(C1+1);     % index of v2*        
           nums3 = nums2+1;            % index of v3*
           nums4 = nums1+1;            % index of v4*
           nume1 = (C1+1)*(cx2-1)+cy2; % index of v1*       
           nume2 = nume1+(C1+1);     % index of v2*        
           nume3 = nume2+1;            % index of v3*
           nume4 = nume1+1;            % index of v4*  
           sp_ii(16*kk-15:16*kk)=[ones(1,8).*(2*kk-1), ones(1,8).*(2*kk)];
           sp_jj(16*kk-15:16*kk)=[2*nums1-1, 2*nums2-1, 2*nums3-1, 2*nums4-1,...
               2*nums1, 2*nums2, 2*nums3, 2*nums4,...
               2*nume1-1, 2*nume2-1, 2*nume3-1, 2*nume4-1,...
               2*nume1, 2*nume2, 2*nume3, 2*nume4] + 2*num_V*(index_j-1);
           sp_ss(16*kk-15:16*kk) = 1/sqrt(length(tmpindex)).*[a.*coeff_mesh_ps; b.*coeff_mesh_ps; a.*coeff_mesh_pe; b.*coeff_mesh_pe];
           
           sp_ci(2*kk-1:2*kk) = [2*kk-1, 2*kk];
           sp_cj(2*kk-1:2*kk) = [js, js];
           sp_cs(2*kk-1:2*kk) = [1, 1];
           
           kk=kk+1;           
       end
       initial_c(js) = 1/sqrt(length(tmpindex))*c;
       js=js+1;
    end
end

% delete redundant zeros in sp_i/sp_ii, sp_j/sp_jj, sp_s/sp_ss, c_Match and c_ones
sp_i = sp_i(1:16*(k-1)); sp_j = sp_j(1:16*(k-1)); sp_s = sp_s(1:16*(k-1));
cMatch = cMatch(1:2*(k-1));
sp_ii = sp_ii(1:16*(kk-1)); sp_jj = sp_jj(1:16*(kk-1)); sp_ss = sp_ss(1:16*(kk-1));
sp_ci = sp_ci(1:2*(kk-1)); sp_cj=sp_cj(1:2*(kk-1)); sp_cs = sp_cs(1:2*(kk-1));
initial_c = initial_c(1:js-1);
sparse_line_ref = sparse(sp_i, sp_j, sp_s, max(sp_i), 2*(num_imgs-1)*num_V);
if isempty(sp_ii)
    c_ones = [];
    sparse_line_free=[];
else
    c_ones = sparse(sp_ci, sp_cj, sp_cs, max(sp_ci), js-1);  
    sparse_line_free = sparse(sp_ii, sp_jj, sp_ss, max(sp_ii), 2*(num_imgs-1)*num_V);
end

% % final sparse matrix of line alignment
% sparse_line_align = [sparse_line_ref, zeros(max(sp_i), size(c_ones,2));
%     sparse_line_free, c_ones];
% cMatch_ones = [cMatch; zeros(size(c_ones, 1),1)];

%% final alignment term (with point and line)
sparse_point_align = [sparse_point_ref, zeros(size(sparse_point_ref,1),size(p_ones,2)), zeros(size(sparse_point_ref,1), size(c_ones,2));
    sparse_point_free, p_ones, zeros(size(sparse_point_free,1), size(c_ones,2))];
sparse_line_align = [sparse_line_ref, zeros(size(sparse_line_ref,1),size(p_ones,2)), zeros(size(sparse_line_ref,1),size(c_ones,2));
    sparse_line_free, zeros(size(sparse_line_free,1),size(p_ones,2)), c_ones];
pMatch_ones = [pMatch; zeros(size(p_ones, 1),1)];
cMatch_ones = [cMatch; zeros(size(c_ones, 1),1)];
initial_pc = [initial_p ; initial_c];

    
end

