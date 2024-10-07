function [sp_cubes] = SpRectangle(sp_label, num_label, img_cube) 
    for idx = 0:num_label-1   % �������г�����
        [row_idx, col_idx] = find(sp_label == idx);	% �ҳ���ǰ�����ص��������ص�����
        row_start = min(row_idx);
        row_end = max(row_idx);
        col_start = min(col_idx);
        col_end = max(col_idx);
        
        % ��ȡ����������������������ݺͳ����ر�ǩ
        temp = img_cube(row_start:row_end, col_start:col_end, :);
        sp_label_cut = sp_label(row_start:row_end, col_start:col_end);
        
        % ���������в����ڵ�ǰ�����ص����ػҶ�ֵ��Ϊ NaN
        [rows, cols, pages] = size(temp);
        temp = reshape(temp, rows * cols, pages);
        loc = find(sp_label_cut(:) ~= idx);
        temp(loc, :) = NaN;
        
        sp_cubes{idx+1} = reshape(temp, rows, cols, pages);
    end
end