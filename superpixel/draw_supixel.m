function draw_supixel(labels,pca3image,fig_num)
[w,h]=size(labels);
dim=size(pca3image,3);
%% normalrite pca3image to (0,255)
image=reshape(double(pca3image),[w*h,dim])';
image=mapminmax(image,0,1)*255;
pca3image=reshape(image',[w h dim]);

%%
% bmap=seg2bmap(labels,h,w);
% bmapOnImg = pca3image;
% grey_img = double(rgb2gray(pca3image));
% idx = find(bmap>0);
% timg = grey_img;
% timg(idx) = 255;
% bmapOnImg(:,:,2) = timg;
% bmapOnImg(:,:,1) = grey_img;
% bmapOnImg(:,:,3) = grey_img;
% h1=figure(1);
% imagesc(bmapOnImg);
% axis('equal','off');
% print(h1,'-depsc','-r300', 'figure\pca_greesu');
%%
h2=figure;
title('Original Image');
imshow(uint8(pca3image));
axis('equal','off');
%title('ground true');
%print(h2,'-depsc','-r300', 'figure\pca_1');
%print(h2,'-djpeg', 'figure\pca_s3');

test=ones(w,h);
for i=1:w-1
    for j=1:h-1
        if labels(i,j)~=labels(i,j+1)
            if test(i,j+1)==1&test(i,j)==1
                test(i,j+1)=0;
            end
        end
        if labels(i,j)~=labels(i+1,j)
            if test(i+1,j)==1&test(i,j)==1
                test(i+1,j)=0;
            end
        end
        if labels(i,j)~=labels(i+1,j+1)
            if test(i+1,j+1)==1&test(i,j)==1
                test(i+1,j+1)=0;
            end
        end
    end
end
%test=uint8(test);
% figure(2);
% imagesc(test);
%bee=imread('data\GT_17etiquetas','pgm');\
k=0;
bee=double(pca3image);
if (dim==3)
    %
    for i=1:size(pca3image,3)
        bee(:,:,i)=bee(:,:,i).*test;
    end
    %
    for i=1:size(pca3image,3)
        newbee=bee(:,:,i);
        %red
        if(i==1)
            newbee(find(test==k))=400;
%             newbee(find(labels==388))=400;
%             newbee(find(labels==387))=0;
        else
            newbee(find(test==k))=0;
%             newbee(find(labels==388))=0;
%             newbee(find(labels==387))=400;
        end
        bee(:,:,i)=newbee;
    end
end
%
if (dim==1)
    bee=bee.*test;
    bee(find(test==k))=400;
end

h3=figure;
title('Superpixel Segmentation Visualization');
imshow(uint8(bee));
axis('equal','off');
% title('super frome pca');
% print(h3,'-djpeg','figure\24148_70_20_super');

%% ground true super map
% if(dim==1)
%     ground_super=label_data;
%     ground_super(find(test==k))=17;
%     h4=figure(4);
%     imagesc(ground_super);
%     axis('equal','off');
%     colormap_rgb = ...
%         [
%         1 1 1;  % white
%         0 0 0;  % black
%         1 0 0;  % red
%         0 1 0;  % green
%         0 0 1;  % blue
%         1 1 0;  % yellow
%         1 0 1;  % carmine
%         0 1 1;  % cyanine
%         2/3 0 1; % sky blue
%         1 0 2/3;
%         0 2/3 1;
%         1 0.5 0; % orange
%         0.5 1 0;
%         0 0 0.5;
%         0.5 0 0; % deep red
%         0 0.5 0;
%         0.5 0.5 0.5; % gray
%         0.1 0.1 0.1; %super line
%         ];
%     rowStart=2;
%     colStart=2;
%     colormap(colormap_rgb);
%     rectangle('Position',[colStart-1,rowStart-1,size(label_data,2)-1,size(label_data,1)-1],'LineWidth',1.5);
% %     print(h4,'-djpeg','figure\pca_1_groundtrue_super');
% end
%%

% label_data(find(label_data==0))=20;
% figure(4);
% imagesc(label_data);
% label_data=uint8(label_data).*test;
% figure(5);
% imagesc(label_data);
end
%}