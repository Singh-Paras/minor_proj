function [BoundInp,BoundLoc,Contourout, Contourin] = getBoundary(rows, cols, mask)
% �˺������������صı߽�
% Parameters��
%   rows:ͼ����
%   cols��ͼ����
%   mask��rows��cols�еĶ�ֵ����
% Return:   BoundInp: ��߽��Ӧ���к��У����е�һ�б�ʾ�У��ڶ��б�ʾ��
%           BoundLoc: �ڱ߽��Ӧ���к��У����е�һ�б�ʾ�У��ڶ��б�ʾ��
%           Contourout: ��߽�
%           Contourin: �ڱ߽�

[row,col] = ndgrid(1:rows,1:cols);
% row = int32(row);
% col = int32(col);
%img = double(ImgData == 0);
img = 1- mask;
% numel(find(img(:)==1));
if(sum(sum(img)) ~= 0)
    % Find out the edge of missing data region.
    se = strel('square',3);
    % ���ͺ�Ϊ1�����򣬱߽����������һ�����ء�
    dilImg = imdilate(img,se);
    % ���ͺ��ͼ���ȥԭͼ�������߽�
    Contourout = dilImg - img;
  
    % ��ʴ��Ϊ1�����򣬱߽����������һ������
    eroImg = imerode(img,se);
    % ��ʴ���ͼ���ȥԭͼ������ڱ߽�
    Contourin = img - eroImg;
    %% ��߽�����
    BoundInp_row_Loc = row .* Contourout;
    BoundInp_col_Loc = col .* Contourout;
    %% ��������
    BoundInp(:,1) = BoundInp_row_Loc(BoundInp_row_Loc ~= 0);
    BoundInp(:,2) = BoundInp_col_Loc(BoundInp_col_Loc ~= 0);
    
    %% �ڱ߽�����
    BoundLoc_row_Loc = row .* Contourin;
    BoundLoc_col_Loc = col .* Contourin;
    %% ��������
    BoundLoc(:,1) = BoundLoc_row_Loc(BoundLoc_row_Loc ~= 0);
    BoundLoc(:,2) = BoundLoc_col_Loc(BoundLoc_col_Loc ~= 0);   
else
    BoundInp = [];
    BoundLoc = [];
    Contourout = zeros(rows,cols);
end
end