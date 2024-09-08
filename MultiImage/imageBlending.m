function [ output ] = imageBlending( cell_Homo, cell_HomoMask, blend_type )

img_n = size(cell_Homo,1);

% composite by average blending
if strcmp(blend_type, 'average')==1
    output = zeros(size(cell_Homo{1}));
    % blend warped images
    sum_ws = zeros(size(output,1), size(output,2));
    for i=1:img_n
        tmp_w = cell_HomoMask{i};
        sum_ws = sum_ws + tmp_w;
        output(:,:,1) = output(:,:,1) + cell_Homo{i}(:,:,1).*tmp_w;
        output(:,:,2) = output(:,:,2) + cell_Homo{i}(:,:,2).*tmp_w;
        output(:,:,3) = output(:,:,3) + cell_Homo{i}(:,:,3).*tmp_w;
    end
    output = output./cat(3, sum_ws, sum_ws, sum_ws);
end

% composite by linear blending
if strcmp(blend_type, 'linear')==1
    for i = 1:img_n
        if i == 1
            out = cell_Homo{1};
            out_mask = cell_HomoMask{1};
            % center of out
            [r, c] = find(cell_HomoMask{1});
            out_center = [mean(r) mean(c)];
        else % blend out and c1out{i}
            % center of c1out{i}
            [r, c] = find(cell_HomoMask{i});
            out_i_center = [mean(r) mean(c)];
            % compute weighted mask
            vec = out_i_center - out_center; % vector from out_center to out_i_center
            intsct_mask = cell_HomoMask{i} & out_mask; % 1 channel

            [r, c] = find(intsct_mask);
            idx = sub2ind(size(cell_HomoMask{i}), r, c);
            out_wmask = zeros(size(cell_HomoMask{i}));
            proj_val = (r - out_center(1))*vec(1) + (c- out_center(2))*vec(2); % inner product
            out_wmask(idx) = (proj_val - (min(proj_val)+(1e-3))) / ...
                             ((max(proj_val)-(1e-3)) - (min(proj_val)+(1e-3))); % weight map (of overlapped area) for c1out{i}, 1 channel
            % blending
            mask1 = out_mask&(out_wmask==0);
            mask2 = out_wmask;
            mask3 = cell_HomoMask{i}&(out_wmask==0);
            mask1 = cat(3, mask1, mask1, mask1);
            mask2 = cat(3, mask2, mask2, mask2);
            mask3 = cat(3, mask3, mask3, mask3);
            out = out.*(mask1+(1-mask2).*(mask2~=0)) + cell_Homo{i}.*(mask2+mask3);
            % update
            out_mask = out_mask | cell_HomoMask{i};
            out_center = out_i_center; % update out_center by assign center of c1out{i}
        end
    end
    output=out;
end

end

