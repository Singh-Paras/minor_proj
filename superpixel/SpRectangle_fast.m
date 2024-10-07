function [sp_cubes] = SpRectangle_fast(sp_label, num_label, img_cube, max_row_length, max_col_length, i)
    pages = size(img_cube, 3);
    sp_cubes = zeros(num_label, max_row_length, max_col_length, pages);
    sp_cubes(:) = NaN;
    for idx = 0:num_label-1   % 遍历所有超像素
        [row_idx, col_idx] = find(sp_label == idx);	% 找出当前超像素的所有像素的索引
        row_start = min(row_idx);
        row_end = max(row_idx);
        col_start = min(col_idx);
        col_end = max(col_idx);
        
        % 提取超像素所属矩形区域的数据和超像素标签
        temp = img_cube(row_start:row_end, col_start:col_end, :);
        sp_label_cut = sp_label(row_start:row_end, col_start:col_end);
        
        % 矩形区域中不属于当前超像素的像素灰度值置为 NaN
        [rows, cols, pages] = size(temp);
        temp = reshape(temp, rows * cols, pages);
        loc = find(sp_label_cut(:) ~= idx);
        temp(loc, :) = NaN;
        
        global band_feature_ERS;
        band_feature_ERS(i, idx+1, 1:rows, 1:cols, :) = reshape(temp, rows, cols, pages);
    end
end