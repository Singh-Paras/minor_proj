function draw_ers_supixel = Untitled2( pca3image,label_ers )
%pca3image    N x M matrics each element of the matrics is between 0~255
%label_ers        ERS�ָ���label����
[w h] = size(label_ers);
test=ones(w,h);%���߽�ͼ(����ʵ��)
for i=1:w-1
    for j=1:h-1
        if label_ers(i,j)~=label_ers(i,j+1)%ͬһ�У����������ڲ�ͬ��
            if test(i,j+1)==1&test(i,j)==1
                test(i,j+1)=0;
            end
        end
        if label_ers(i,j)~=label_ers(i+1,j)%ͬһ�У����������ڲ�ͬ��
            if test(i+1,j)==1&test(i,j)==1
                test(i+1,j)=0;
            end
        end
        if label_ers(i,j)~=label_ers(i+1,j+1)
            if test(i+1,j+1)==1&test(i,j)==1
                test(i+1,j+1)=0;
            end
        end
    end
end

k=0;
bee=double(pca3image);
%
bee=bee.*test;
bee(find(test==k))=400;


h3=figure;
imshow(uint8(bee));
axis('equal','off');
end

