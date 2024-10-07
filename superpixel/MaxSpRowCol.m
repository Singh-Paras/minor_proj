function [max_row_length, max_col_length] = MaxSpRowCol(sp_label, num_label) 
    max_row_length = 0;
    max_col_length = 0;
    for idx = 0:num_label-1   % 遍历所有超像素
        [row_idx, col_idx] = find(sp_label == idx);	% 找出当前超像素的所有像素的索引
        row_start = min(row_idx);
        row_end = max(row_idx);
        col_start = min(col_idx);
        col_end = max(col_idx);
        
        if row_end - row_start + 1 > max_row_length
            max_row_length = row_end - row_start + 1;
        end
        if col_end - col_start + 1 > max_col_length
            max_col_length = col_end - col_start + 1;
        end
    end
end