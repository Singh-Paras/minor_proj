clc; clear; close all; warning('off');


%% 加载高光谱数据
load Salinas_corrected.mat;
load Salinas_pca_3.mat;

cube_data = double(salinas_corrected);
[rows, cols, pages] = size(cube_data);
cube_data_flat = reshape(cube_data, rows * cols, pages)';	% 三维变二维
cube_data_nm_flat = mapminmax(cube_data_flat, 0, 255) ;  % 数据归一化0~255
cube_data_nm = reshape(cube_data_nm_flat', rows, cols, pages);

T = 3.7; % Salinas  
K = [1,2,3];
Nk = zeros(1,length(K));
sp_label_SLIC = zeros(rows, cols, length(Nk));
num_label_SLIC = zeros(1, length(Nk));


%% 超像素分割 SLIC
for i=1:length(K)
    Nk(i) = floor(K(i)*(rows*cols)/(T*T));
    [sp_label_SLIC(:,:,i), num_label_SLIC(i)]= slicomex(Salinas_pca_3, Nk(i));	% SLIC 超像素分割
    num_label_SLIC(i) = max(max(sp_label_SLIC(:,:,i)))+1;
%     draw_supixel(sp_label_SLIC(:,:,i), cube_data(:, :, 1:3), 1);	% 超像素分割可视化
   
    %% 计算超像素之间的距离
 
    band_feature_SLIC{i} = SpRectangle(sp_label_SLIC(:,:,i), num_label_SLIC(i) , cube_data_nm);
    sim_KL_SLIC{i} = zeros(pages,pages,num_label_SLIC(i));
    sim_L2_SLIC{i} = zeros(pages,pages,num_label_SLIC(i));
    sim_L21_SLIC{i} = zeros(pages,pages);
    sim_KL1_SLIC{i} = zeros(pages,pages);

    for t=1:num_label_SLIC(i)
        for j = 1:pages-1
            Vj = band_feature_SLIC{i}{t}(:,:,j);
            Vj(isnan(Vj)) = [];
            for k = j+1:pages
                Vk =band_feature_SLIC{i}{t}(:,:,k);
                Vk(isnan(Vk)) = [];
                sim_KL_SLIC{i}(j,k,t) = KL(Vj(:)'+eps,Vk(:)'+eps);
                sim_L2_SLIC{i}(j,k,t) = norm((Vj-Vk),2)/pages;
                sim_KL_SLIC{i}(k,j,t) = sim_KL_SLIC{i}(j,k,t);
                sim_L2_SLIC{i}(k,j,t) = sim_L2_SLIC{i}(j,k,t);
            end
        end
    end
    
    for j = 1:pages-1
        for k = j+1:pages
            temp2 = sim_L2_SLIC{i}(j,k,:);
            sim_L21_SLIC{i}(j,k) = norm(temp2(:),1);
            sim_L21_SLIC{i}(k,j) =  sim_L21_SLIC{i}(j,k);
        end
    end
    
    for j = 1:pages-1
        for k = j+1:pages
            temp2 = sim_KL_SLIC{i}(j,k,:);
            sim_KL1_SLIC{i}(j,k) = norm(temp2(:),1);
            sim_KL1_SLIC{i}(k,j) =  sim_KL1_SLIC{i}(j,k);
        end
    end
  
%     figure;
%     imagesc(sim_L21_SLIC{i});
%     figure;
%     imagesc(sim_KL1_SLIC{i}); 
end

save('Salinas_sim_L21_SLIC','sim_L21_SLIC');
save('Salinas_sim_KL1_SLIC','sim_KL1_SLIC');