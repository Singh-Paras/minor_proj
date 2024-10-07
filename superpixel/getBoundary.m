function [BoundInp,BoundLoc,Contourout, Contourin] = getBoundary(rows, cols, mask)
% 此函数用来求超像素的边界
% Parameters：
%   rows:图像行
%   cols：图像列
%   mask：rows行cols列的二值矩阵
% Return:   BoundInp: 外边界对应的行和列，其中第一列表示行，第二列表示列
%           BoundLoc: 内边界对应的行和列，其中第一列表示行，第二列表示列
%           Contourout: 外边界
%           Contourin: 内边界

[row,col] = ndgrid(1:rows,1:cols);
% row = int32(row);
% col = int32(col);
%img = double(ImgData == 0);
img = 1- mask;
% numel(find(img(:)==1));
if(sum(sum(img)) ~= 0)
    % Find out the edge of missing data region.
    se = strel('square',3);
    % 膨胀后为1的区域，边界会向外扩张一个像素。
    dilImg = imdilate(img,se);
    % 膨胀后的图像减去原图像就是外边界
    Contourout = dilImg - img;
  
    % 腐蚀后为1的区域，边界会向内收缩一个像素
    eroImg = imerode(img,se);
    % 腐蚀后的图像减去原图像就是内边界
    Contourin = img - eroImg;
    %% 外边界坐标
    BoundInp_row_Loc = row .* Contourout;
    BoundInp_col_Loc = col .* Contourout;
    %% 拉成向量
    BoundInp(:,1) = BoundInp_row_Loc(BoundInp_row_Loc ~= 0);
    BoundInp(:,2) = BoundInp_col_Loc(BoundInp_col_Loc ~= 0);
    
    %% 内边界坐标
    BoundLoc_row_Loc = row .* Contourin;
    BoundLoc_col_Loc = col .* Contourin;
    %% 拉成向量
    BoundLoc(:,1) = BoundLoc_row_Loc(BoundLoc_row_Loc ~= 0);
    BoundLoc(:,2) = BoundLoc_col_Loc(BoundLoc_col_Loc ~= 0);   
else
    BoundInp = [];
    BoundLoc = [];
    Contourout = zeros(rows,cols);
end
end