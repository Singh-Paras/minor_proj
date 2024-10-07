function [neighbor_supixel, used_time] = supixelNeighbor(supixel_num, label_ers, rows, cols)
neighbor_supixel = cell(1, supixel_num);
tstart = tic;
for index = 1:supixel_num
    mask = label_ers + 1;
    mask(mask == index) = 0;
    mask(mask ~= 0) = 1;
    [~,~,Contourout, ~] = getBoundary(rows, cols, mask);
    Contourout_label = (Contourout==1);
    near_supixel_temp = unique(label_ers(Contourout_label));
    neighbor_supixel{1, index} = near_supixel_temp;
    near_supixel_temp = [];
    mask = [];
end
neighbor_supixel = neighbor_supixel';
used_time = toc(tstart);
end