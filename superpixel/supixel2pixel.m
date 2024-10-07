function [result] = supixel2pixel(supixel_result, label_ers, tr_index, un_index)
result = zeros(1, length(label_ers));
index = [tr_index un_index];

for i = 1:length(index)
    pos = (label_ers == (index(i)));
    result(pos) = supixel_result(i);
end
end