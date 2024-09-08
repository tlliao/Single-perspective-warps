function [ baH, ref_img, points_xk, points_predict, lines_xk, lines_predict] = bundleAdjustmentRefofLine(data_orig, inliers, data_lines, h, imgs_pairs)
% given pairwise inliers, original matches and valid image pairs, use BA
% algorithm to calculate homographies r.t. reference
%fprintf(' Performing Bundle Adjustment...');tic;

num_imgs = max(imgs_pairs(:,2));
valid_pairs = find(imgs_pairs(:,3));

% extract reference image and stitching order.
%---------------------------------------------
% Build image connection graph. This graph will be useful for 
% generating the stitching order and determining the reference image.
%fprintf('>Building image connection graph...');tic;
imgs_conn_graph = zeros(num_imgs,num_imgs);
for i=1:size(valid_pairs,1)
    % The weight of each arc in the graph is the number of (pairwise) keypoint matches.
    imgs_conn_graph(imgs_pairs(valid_pairs(i),1),imgs_pairs(valid_pairs(i),2)) = length(inliers{valid_pairs(i)});
    imgs_conn_graph(imgs_pairs(valid_pairs(i),2),imgs_pairs(valid_pairs(i),1)) = length(inliers{valid_pairs(i)});
end
%fprintf('done (%fs)\n',toc);

% choose reference image (based on image graph).
%fprintf('>Obtaining reference image and stitching order...');tic;
[~, ref_img] = max(sum(imgs_conn_graph,2));
if abs((ref_img-(1+num_imgs)/2))>0.5; ref_img = floor((1+num_imgs)/2); end


% find stitching order.
stitching_order = {ref_img};
ordered_imgs = ref_img; % the first ordered image is the reference image.
unordered_imgs = setdiff(1:num_imgs,ordered_imgs); % unordered images.
graph_paths = cell(num_imgs,1); % contains the path from one image to the reference image.
graph_paths{ref_img} = ref_img;
while ~isempty(unordered_imgs)    
    % Select the next image (using the weighted image connection graph).
    edge_weight = 0;
    for i=unordered_imgs
        for j=ordered_imgs
            if imgs_conn_graph(i,j) > edge_weight
                edge_weight = imgs_conn_graph(i,j);
                next_img   = i;
                middle_img = j;
            end
        end
    end
    
    if edge_weight == 0 % There was no next image, we probably
                        % have a disconnected panorama.
        fprintf('\n  ');
        warning('Found a disconnected panorama...');
        break;
    end

    % Find the path from the next image to the reference image.
    graph_paths{next_img} = [graph_paths{middle_img} next_img];
    stitching_order = [stitching_order; {graph_paths{next_img}}];
    
    % Update variables.            
    ordered_imgs = [ordered_imgs next_img];
    unordered_imgs = setdiff(1:num_imgs,ordered_imgs);    
end
%fprintf('done (%fs)\n',toc);

% Order data according to the stitching order (target image is always on the left).
%fprintf('> Ordering data according to stitching order...');tic;
stitch_pairs = []; % Will save the (valid) image pairs that we need for generating 
                   % the panorama (following the stitching order).
for i=1:size(stitching_order,1)
    for j=1:size(stitching_order{i},2)-1        
        for k=1:size(valid_pairs,1)
            if stitching_order{i}(j) == imgs_pairs(valid_pairs(k),2) && stitching_order{i}(j+1) == imgs_pairs(valid_pairs(k),1)                
                aux = imgs_pairs(valid_pairs(k),1);
                imgs_pairs(valid_pairs(k),1) = imgs_pairs(valid_pairs(k),2);
                imgs_pairs(valid_pairs(k),2) = aux;

                aux = data_orig{valid_pairs(k)}(1:3,:);
                data_orig{valid_pairs(k)}(1:3,:) = data_orig{valid_pairs(k)}(4:6,:);
                data_orig{valid_pairs(k)}(4:6,:) = aux;
                
                aux = data_lines{valid_pairs(k)}(:,1:4);
                data_lines{valid_pairs(k)}(:,1:4) = data_lines{valid_pairs(k)}(:,5:8);
                data_lines{valid_pairs(k)}(:,5:8) = aux;
                
                h{valid_pairs(k)} = reshape(reshape(h{valid_pairs(k)},3,3)\eye(3),9,1);
                
                stitch_pairs = [stitch_pairs; valid_pairs(k)];
                break;
            elseif stitching_order{i}(j) == imgs_pairs(valid_pairs(k),1) && stitching_order{i}(j+1) == imgs_pairs(valid_pairs(k),2)
                stitch_pairs = [stitch_pairs; valid_pairs(k)];
                break;
            end
        end
    end
end
stitch_pairs = unique(stitch_pairs);
%fprintf('done (%fs)\n',toc);

%------------------------------------------------------
% Image stitching with pairwise (chained) Homographies.
%------------------------------------------------------
%fprintf('Image stitching with pairwise Homographies\n');
% Chaining the homographies H (obtained by means of
% pairwise DLT). These chained H's will be the initial 
% estimates for BA.
H = cell(size(stitching_order,1),1);
for i=1:size(stitching_order,1)
    H{i} = eye(3);
    for j=1:size(stitching_order{i},2)-1
        for k=1:size(stitch_pairs,1)
            if stitching_order{i}(j) == imgs_pairs(stitch_pairs(k),1) && stitching_order{i}(j+1) == imgs_pairs(stitch_pairs(k),2)
                H{i} = reshape(h{stitch_pairs(k)},3,3) * H{i};
                break;
            end
        end
    end
end

% sort the paramaters following the image order.
lm_params0  = [];
for i=1:num_imgs
    if i==ref_img; continue; end  % there is no need to consider reference
    for j=1:size(stitching_order,1)
        if i == stitching_order{j}(end)
            lm_params0 = [lm_params0; H{j}(:)];
            break;
        end
    end
end

%-----------------------------------------
% Image stitching using Bundle Adjustment.
%-----------------------------------------
% Obtaining size of canvas.
%fprintf('Image stitching with Bundle Adjustment\n');
% Obtain xk's on images
%fprintf('> Getting xk''s on images...');tic;
all_kpts_per_img = cell(num_imgs,1);
all_klines_per_img = cell(num_imgs,1);
for i=1:num_imgs
    for j=1:size(stitch_pairs,1)
        if imgs_pairs(stitch_pairs(j),1) == i
            all_kpts_per_img{i} = unique([all_kpts_per_img{i} data_orig{stitch_pairs(j)}(1:2,inliers{stitch_pairs(j)})]','rows')';
            all_klines_per_img{i} = unique([all_klines_per_img{i}; data_lines{stitch_pairs(j)}(:,1:4)],'rows');
        elseif imgs_pairs(stitch_pairs(j),2) == i
            all_kpts_per_img{i} = unique([all_kpts_per_img{i} data_orig{stitch_pairs(j)}(4:5,inliers{stitch_pairs(j)})]','rows')';
            all_klines_per_img{i} = unique([all_klines_per_img{i}; data_lines{stitch_pairs(j)}(:,5:8)],'rows');
        end
    end
end
%fprintf('done (%fs)\n',toc);

% Obtain xi's on canvas
%fprintf('  Getting xi''s on canvas...');tic;
xk = []; 
for i=1:num_imgs
    for j=1:size(all_kpts_per_img{i},2)
        aux = NaN(2*(num_imgs),1);
        aux(((i-1)*2)+1:i*2) = all_kpts_per_img{i}(:,j);
        for k=1:size(stitch_pairs,1)
            if imgs_pairs(stitch_pairs(k),1) == i
                idx = getLocation(data_orig{stitch_pairs(k)}(1:2,inliers{stitch_pairs(k)})',all_kpts_per_img{i}(:,j)');
                if idx ~= 0
                    aux(((imgs_pairs(stitch_pairs(k),2)-1)*2)+1:imgs_pairs(stitch_pairs(k),2)*2) = data_orig{stitch_pairs(k)}(4:5,inliers{stitch_pairs(k)}(idx));
                end
            elseif imgs_pairs(stitch_pairs(k),2) == i
                idx = getLocation(data_orig{stitch_pairs(k)}(4:5,inliers{stitch_pairs(k)})',all_kpts_per_img{i}(:,j)');
                if idx ~= 0
                    aux(((imgs_pairs(stitch_pairs(k),1)-1)*2)+1:imgs_pairs(stitch_pairs(k),1)*2) = data_orig{stitch_pairs(k)}(1:2,inliers{stitch_pairs(k)}(idx));
                end
            end
        end
        if ~isempty(xk)
            found = 0;
            for k=1:num_imgs
                idx = getLocation(xk((k-1)*2+1:k*2,:)',aux((k-1)*2+1:k*2)');
                if idx ~= 0
                    found = 1;
                    for l=1:num_imgs
                        if ~isnan(aux((l-1)*2+1))
                            xk((l-1)*2+1:l*2,idx) = aux((l-1)*2+1:l*2);
                        end
                    end
                    break;
                end
            end
            if found == 0
                xk = [xk aux];
            end
        else
            xk = [xk aux];
        end
    end
end

lm_xi = zeros(3,size(xk,2));
xik_idx = cell(num_imgs,1);
num_xik = zeros(num_imgs,1);
for i=1:num_imgs
    for j=1:size(stitching_order,1)
        if i == stitching_order{j}(end)
            xik_idx{i} = find(~isnan(xk((i-1)*2+1,:)));
            aux = H{j} \ [xk((i-1)*2+1:i*2,xik_idx{i});ones(1,size(xik_idx{i},2))];
            lm_xi(1:3,xik_idx{i}) = lm_xi(1:3,xik_idx{i}) +...
                [aux(1,:)./aux(3,:); aux(2,:)./aux(3,:); ones(1,size(xik_idx{i},2))];
            break;
        end
    end
    num_xik(i) = size(xik_idx{i},2);
end

lm_xi(1,:) = lm_xi(1,:) ./ lm_xi(3,:);
lm_xi(2,:) = lm_xi(2,:) ./ lm_xi(3,:);
%fprintf('done (%fs)\n',toc);

% Obtain lines_xi's on canvas
%fprintf('  Getting lines_xi''s on canvas...');tic;
lines_xk = []; 
for i=1:num_imgs
    for j=1:size(all_klines_per_img{i},1)
        aux = NaN(4*(num_imgs),1);
        aux(4*i-3:4*i) = all_klines_per_img{i}(j,:)';
        for k=1:size(stitch_pairs,1)
            if imgs_pairs(stitch_pairs(k),1) == i
                idx = linesFindLocation(data_lines{stitch_pairs(k)}(:,1:4), all_klines_per_img{i}(j,:));
                if idx ~= 0
                    aux(4*imgs_pairs(stitch_pairs(k),2)-3:4*imgs_pairs(stitch_pairs(k),2)) = data_lines{stitch_pairs(k)}(idx,5:8)';
                end
            elseif imgs_pairs(stitch_pairs(k),2) == i
                idx = linesFindLocation(data_lines{stitch_pairs(k)}(:,5:8), all_klines_per_img{i}(j,:));
                if idx ~= 0
                    aux(4*imgs_pairs(stitch_pairs(k),1)-3:4*imgs_pairs(stitch_pairs(k),1)) = data_lines{stitch_pairs(k)}(idx,1:4)';
                end
            end
        end
        if ~isempty(lines_xk)
            found = 0;
            for k=1:num_imgs
                idx = linesFindLocation(lines_xk(4*k-3:4*k,:)', aux(4*k-3:4*k)'); 
                if idx ~= 0
                    found = 1;
                    for l=1:num_imgs
                        if ~isnan(aux(4*l-3))
                            lines_xk(4*l-3:4*l,idx) = aux(4*l-3:4*l);
                        end
                    end
                    break;
                end
            end
            if found == 0
                lines_xk = [lines_xk aux];
            end
        else
            lines_xk = [lines_xk aux];
        end
    end
end

lines_lm_xi = zeros(6,size(lines_xk,2));
lines_xik_idx = cell(num_imgs,1);
lines_num_xik = zeros(num_imgs,1);
for i=1:num_imgs
    for j=1:size(stitching_order,1)
        if i == stitching_order{j}(end)
            lines_xik_idx{i} = find(~isnan(lines_xk(4*i-3,:)));
            aux1 = H{j} \ [lines_xk(4*i-3:4*i-2,lines_xik_idx{i});ones(1,size(lines_xik_idx{i},2))];
            aux2 = H{j} \ [lines_xk(4*i-1:4*i,lines_xik_idx{i});ones(1,size(lines_xik_idx{i},2))];
            lines_lm_xi(1:3,lines_xik_idx{i}) = lines_lm_xi(1:3,lines_xik_idx{i}) +...
                [aux1(1,:)./aux1(3,:); aux1(2,:)./aux1(3,:); ones(1,size(lines_xik_idx{i},2))];
            lines_lm_xi(4:6,lines_xik_idx{i}) = lines_lm_xi(4:6,lines_xik_idx{i}) +...
                [aux2(1,:)./aux2(3,:); aux2(2,:)./aux2(3,:); ones(1,size(lines_xik_idx{i},2))];
            break;
        end
    end
    lines_num_xik(i) = size(lines_xik_idx{i},2);
end

lines_lm_xi(1,:) = lines_lm_xi(1,:) ./ lines_lm_xi(3,:);
lines_lm_xi(2,:) = lines_lm_xi(2,:) ./ lines_lm_xi(3,:);
lines_lm_xi(4,:) = lines_lm_xi(4,:) ./ lines_lm_xi(6,:);
lines_lm_xi(5,:) = lines_lm_xi(5,:) ./ lines_lm_xi(6,:);
%fprintf('done (%fs)\n',toc);

%fprintf('  Performing Bundle Adjustment (Levenberg-Marquardt)...');
% Set of parameters to be optimised (H1,...,Hk,x1,...,xN)
lm_params0   = [lm_params0; reshape(lm_xi(1:2,:),size(lm_xi(1:2,:),2)*2,1); reshape(lines_lm_xi([1,2,4,5],:),size(lines_lm_xi,2)*4,1)];

% Normalizer (lm_norm contains |M(i)| per point xi)
lm_norm = lm_xi(3,:)';  lines_lm_norm = lines_lm_xi(6,:)';

% Get the set of observations (xik)
lm_xi_x = [];
lm_xi_y = [];
for i=1:num_imgs
    lm_xi_x = [lm_xi_x xk((i-1)*2+1, xik_idx{i})];
    lm_xi_y = [lm_xi_y xk((i-1)*2+2, xik_idx{i})];    
end
lm_xk = [lm_xi_x lm_xi_y]';

% generate ax+by+c=0 for each lines in lines_xk
lines_abc = [];
for i=1:num_imgs
    aux = lines_xk(4*i-3:4*i,lines_xik_idx{i});
    sqrt_d = ((aux(4,:)-aux(2,:)).^2 + (aux(1,:)-aux(3,:)).^2).^(1/2);
    aux_a = (aux(4,:)-aux(2,:))./sqrt_d;
    aux_b = (aux(1,:)-aux(3,:))./sqrt_d;
    aux_c = (aux(3,:).*aux(2,:)-aux(1,:).*aux(4,:))./sqrt_d;
    lines_abc = [lines_abc; aux_a'; aux_b'; aux_c'];   
end
  
% Perform Bundle Adjustment (sparse Levenberg-Marquardt).
lm_params = ceresPointLineofLengthRef(double(lm_params0), xik_idx, lm_norm, num_imgs, ref_img, lm_xk, lines_xik_idx, lines_lm_norm, 50.*lines_abc);
% fprintf('done (%fs)\n',toc);

baH = cell(num_imgs,1);  k=1;
for i=1:size(stitching_order,1)
    if i==ref_img; baH{i}=eye(3); continue; end
    baH{i} = reshape(lm_params(((k-1)*9)+1:k*9),3,3);
    k=k+1;
end

%fprintf('done (%fs)\n',toc); 

% make the features used in mesh deformation are same with bundle
% adjustment  by ltl
points_xk = xk;
points_predict = reshape(lm_params((num_imgs-1)*9+1:(num_imgs-1)*9+2*length(lm_norm)), 2, length(lm_norm));
points_predict(:, xik_idx{ref_img}) = points_xk(2*ref_img-1:2*ref_img, xik_idx{ref_img});

% calculate the predicted lines functions  ax+by+c=0  (c unknown)
lines_predict = reshape(lm_params((num_imgs-1)*9+2*length(lm_norm)+1:end),4,length(lines_lm_norm));
lines_predict(:, lines_xik_idx{ref_img}) = lines_xk(4*ref_img-3:4*ref_img, lines_xik_idx{ref_img});


end

