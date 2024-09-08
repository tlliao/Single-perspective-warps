function [ sparse_v ] = energyLineV(cell_imgs, C1, C2, cell_v, nor_vecs)
% for all v sample points in target images, preserve the slope and
% equidistant of all lines
num_imgs = size(cell_imgs,1); % number of total images (include reference)
for i=1:size(cell_v,1)
    if isempty(cell_v{i})
        ref_n=i; break;  % the index of reference
    end
end 
%% generating energy v-line term (slope-invariant and equidistant)
Mesh_ps = zeros(4,2); Mesh_pm = zeros(4,2); Mesh_pe = zeros(4,2);
num_V = (C1+1)*(C2+1);  % number of control vertices
k=1; kk=1; row_sp=0; row_spp=0;
for i=[1:ref_n-1,ref_n+1:num_imgs]
    row_sp = row_sp + sum(cell_v{i}(1:2:end-1,end)-1);
    row_spp = row_spp + 2*sum(max(0,cell_v{i}(1:2:end-1,end)-2));
end

%% slope-invariant
sp_i = zeros(16*row_sp,1); sp_j=sp_i; sp_s=sp_i; % three indices of sparse function
for vi=[1:ref_n-1,ref_n+1:num_imgs]
    V_index = (vi<ref_n)*vi+(vi>ref_n)*(vi-1);
    lines_v = cell_v{vi};  nor_vec = nor_vecs(vi,:);
    [M, N, ~] = size(cell_imgs{vi});
    X_col = linspace(1, N, C2+1); % column index of cells
    Y_row = linspace(1, M, C1+1); % row index of cells
    x_dis = X_col(2)-X_col(1);  % the width of scale-cell
    y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell
    for i=1:2:size(lines_v,1)-1
        num_s = lines_v(i,end);  % number of sample points in this segment
        if num_s<=1; continue; end % if sample points less than 2, continue
        for j=1:num_s-1
            lps = [lines_v(i,j),  lines_v(i+1,j)];
            lpe = [lines_v(i,j+1),lines_v(i+1,j+1)];
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
            sp_s(16*k-15:16*k) = [-nor_vec(1).*coeff_mesh_ps; nor_vec(1).*coeff_mesh_pe; -nor_vec(2).*coeff_mesh_ps; nor_vec(2).*coeff_mesh_pe];
            k = k + 1;
        end
    end
end

%% equidistant
sp_ii = zeros(12*row_spp,1); sp_jj=sp_ii; sp_ss=sp_ii; % three indices of sparse function
for vi= [1:ref_n-1,ref_n+1:num_imgs]
    V_index = (vi<ref_n)*vi+(vi>ref_n)*(vi-1);
    lines_v = cell_v{vi}; 
    [M, N, ~] = size(cell_imgs{vi});
    X_col = linspace(1, N, C2+1); % column index of cells
    Y_row = linspace(1, M, C1+1); % row index of cells
    x_dis = X_col(2)-X_col(1);  % the width of scale-cell
    y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell
    for i=1:2:size(lines_v,1)-1
        num_s = lines_v(i,end);  % number of sample points in this segment
        if num_s<=2; continue; end % if sample points less than 2, continue
        for j=2:num_s-1
            lps = [lines_v(i,j-1), lines_v(i+1,j-1)];
            lpm = [lines_v(i,j),   lines_v(i+1,j)];
            lpe = [lines_v(i,j+1), lines_v(i+1,j+1)];
            pxs = min( find(lps(1)-X_col<x_dis & lps(1)-X_col>=0, 1), C2); % the x index of p's position
            pys = min( find(lps(2)-Y_row<y_dis & lps(2)-Y_row>=0, 1), C1); % the y index of p's position
            pxm = min( find(lpm(1)-X_col<x_dis & lpm(1)-X_col>=0, 1), C2); % the x index of p's position
            pym = min( find(lpm(2)-Y_row<y_dis & lpm(2)-Y_row>=0, 1), C1); % the y index of p's position
            pxe = min( find(lpe(1)-X_col<x_dis & lpe(1)-X_col>=0, 1), C2); % the x index of p's position
            pye = min( find(lpe(2)-Y_row<y_dis & lpe(2)-Y_row>=0, 1), C1); % the y index of p's position

            nums1 = (C1+1)*(pxs-1) + pys; % index of v1*   
            nums2 = nums1 + C1+1;
            nums3 = nums2 + 1;
            nums4 = nums1 + 1;
            numm1 = (C1+1)*(pxm-1) + pym;
            numm2 = numm1 + C1+1;
            numm3 = numm2 + 1;
            numm4 = numm1 + 1;
            nume1 = (C1+1)*(pxe-1) + pye;
            nume2 = nume1 + C1+1;
            nume3 = nume2 + 1;
            nume4 = nume1 + 1;        

            Mesh_ps(1:4,:) = [X_col(pxs), Y_row(pys);     % v1
                                X_col(pxs+1), Y_row(pys);   % v2
                                X_col(pxs+1), Y_row(pys+1); % v3
                                X_col(pxs), Y_row(pys+1)];   % v4
            Mesh_pm(1:4,:) = [X_col(pxm), Y_row(pym);     % v1
                                X_col(pxm+1), Y_row(pym);   % v2
                                X_col(pxm+1), Y_row(pym+1); % v3
                                X_col(pxm), Y_row(pym+1)];   % v4      
            Mesh_pe(1:4,:) = [X_col(pxe), Y_row(pye);     % v1
                                X_col(pxe+1), Y_row(pye);   % v2
                                X_col(pxe+1), Y_row(pye+1); % v3
                                X_col(pxe), Y_row(pye+1)];   % v4                        
            coeff_mesh_ps = meshGridAlign(Mesh_ps, lps);
            coeff_mesh_pm = meshGridAlign(Mesh_pm, lpm);
            coeff_mesh_pe = meshGridAlign(Mesh_pe, lpe);
            sp_ii(24*kk-23:24*kk) = [(2*kk-1).*ones(1,12), 2*kk.*ones(1,12)];
            sp_jj(24*kk-23:24*kk) = [2*nums1-1, 2*nums2-1, 2*nums3-1, 2*nums4-1,...
                                  2*numm1-1, 2*numm2-1, 2*numm3-1, 2*numm4-1,...
                                  2*nume1-1, 2*nume2-1, 2*nume3-1, 2*nume4-1,...
                                  2*nums1, 2*nums2, 2*nums3, 2*nums4,...
                                  2*numm1, 2*numm2, 2*numm3, 2*numm4,...
                                  2*nume1, 2*nume2, 2*nume3, 2*nume4] + 2*num_V*(V_index-1);
            sp_ss(24*kk-23:24*kk) = [coeff_mesh_ps; -2.*coeff_mesh_pm; coeff_mesh_pe;...
                                    coeff_mesh_ps; -2.*coeff_mesh_pm; coeff_mesh_pe];
            kk = kk + 1;
        end
    end
end

sparse_v = sparse([sp_i; sp_ii], [sp_j; sp_jj], [sp_s; sp_ss], row_sp + row_spp, 2*(num_imgs-1)*num_V);

end

