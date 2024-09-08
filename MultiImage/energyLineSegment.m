function [ sparse_line ] = energyLineSegment( cell_imgs, cell_lines, cell_slope, baH, C1, C2 )
% given detected line segments, calculate a sparse matrix to preserving
% line struture in image
num_imgs = size(cell_imgs,1); % number of total images (include reference)
for i=1:size(cell_lines,1)
    if isempty(cell_lines{i})
        ref_n=i; break; % the index of reference
    end
end 
num_V = (C1+1)*(C2+1);  % number of control vertices
Mesh_ps = zeros(4,2); Mesh_pe = zeros(4,2);
row_sp=0; k=1;
for i=[1:ref_n-1, ref_n+1:num_imgs]
    row_sp = row_sp + sum(cell_lines{i}(1:2:end-1,end)-1);
end

%% preserving all the detected lines to be lines after transformation
sp_i = zeros(16*row_sp,1); sp_j=sp_i; sp_s=sp_i; % the three indices the sparse function needs
for vi=[1:ref_n-1, ref_n+1:num_imgs]
    V_index = (vi<ref_n)*vi+(vi>ref_n)*(vi-1);
    init_H = baH{ref_n}/baH{vi};   init_H = init_H./init_H(end);
    lines = cell_lines{vi}; slope = cell_slope{vi}; [M, N, ~] = size(cell_imgs{vi});
    X_col = linspace(1, N, C2+1); % column index of cells
    Y_row = linspace(1, M, C1+1); % row index of cells
    x_dis = X_col(2)-X_col(1);  % the width of scale-cell
    y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell
    for i=1:2:size(lines, 1)-1
        num_s = lines(i,end);  % number of sample points in this segment
        if num_s<=1; continue; end
        k_xy = calcSlope(init_H, slope(i), [lines(i,1); lines(i+1,1)]);  % line's slope after transformation
        if isinf(abs(k_xy));  nor_vec=[1,0]; else;  nor_vec = [k_xy, -1];  end
        nor_vec = nor_vec./norm(nor_vec);  % normal vector of warped lines
        for j=1:num_s-1
            lps = [lines(i,j),     lines(i+1,j)];
            lpe = [lines(i,j+1),   lines(i+1,j+1)];
            pxs = min( find(lps(1)-X_col<x_dis & lps(1)-X_col>=0, 1), C2); % the x index of p's position
            pys = min( find(lps(2)-Y_row<y_dis & lps(2)-Y_row>=0, 1), C1); % the y index of p's position
            pxe = min( find(lpe(1)-X_col<x_dis & lpe(1)-X_col>=0, 1), C2); % the x index of p's position
            pye = min( find(lpe(2)-Y_row<y_dis & lpe(2)-Y_row>=0, 1), C1); % the y index of p's position

            nums1 = (C1+1)*(pxs-1) + pys; % index of v1*   
            nums2 = nums1 + C1+1;
            nums3 = nums2 + 1;
            nums4 = nums1 + 1;
            nume1 = (C1+1)*(pxe-1) + pye;
            nume2 = nume1 + C1+1;
            nume3 = nume2 + 1;
            nume4 = nume1 + 1;
            Mesh_ps(1:4,:) = [X_col(pxs), Y_row(pys);     % v1
                                X_col(pxs+1), Y_row(pys);   % v2
                                X_col(pxs+1), Y_row(pys+1); % v3
                                X_col(pxs), Y_row(pys+1)];   % v4
            Mesh_pe(1:4,:) = [X_col(pxe), Y_row(pye);     % v1
                                X_col(pxe+1), Y_row(pye);   % v2
                                X_col(pxe+1), Y_row(pye+1); % v3
                                X_col(pxe), Y_row(pye+1)];   % v4                    
            coeff_mesh_ps = meshGridAlign(Mesh_ps, lps);
            coeff_mesh_pe = meshGridAlign(Mesh_pe, lpe);
            sp_i(16*k-15:16*k) = k.*ones(1,16);
            sp_j(16*k-15:16*k) = [2*nums1-1, 2*nums2-1, 2*nums3-1, 2*nums4-1,...
                                  2*nume1-1, 2*nume2-1, 2*nume3-1, 2*nume4-1,...
                                  2*nums1, 2*nums2, 2*nums3, 2*nums4,...
                                  2*nume1, 2*nume2, 2*nume3, 2*nume4] + 2*num_V*(V_index-1);
            sp_s(16*k-15:16*k) = [-nor_vec(1).*coeff_mesh_ps; nor_vec(1).*coeff_mesh_pe;
                                  -nor_vec(2).*coeff_mesh_ps; nor_vec(2).*coeff_mesh_pe];
            k = k + 1;
        end
    end
end

sparse_line = sparse(sp_i, sp_j, sp_s, row_sp, 2*num_V*(num_imgs-1));

end
