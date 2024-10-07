clc; clear; close all; warning('off');


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
sp_label_ERS = zeros(rows, cols, length(Nk));
num_label_ERS = zeros(1, length(Nk));
pwd
dir
which mymex_ers

%% 超像素分割 ERS

for i=1:length(K)
    Nk(i) = floor(K(i)*(rows*cols)/(T*T));
    sp_label_ERS(:,:,i) = mymex_ers(double(Salinas_pca_3), Nk(i), 0.05, 5);	% ERS 超像素分割 0.05 5
    num_label_ERS(i) = max(max(sp_label_ERS(:,:,i)))+1;
%     draw_supixel(sp_label_ERS(:,:,i), cube_data(:, :, 1:3), 1);	% 超像素分割可视化
   
    %% 计算超像素之间的距离
 
    band_feature_ERS{i} = SpRectangle(sp_label_ERS(:,:,i), num_label_ERS(i) , cube_data_nm);
    sim_KL_ERS{i} = zeros(pages,pages,num_label_ERS(i));
    sim_L2_ERS{i} = zeros(pages,pages,num_label_ERS(i));
    sim_L21_ERS{i} = zeros(pages,pages);
    sim_KL1_ERS{i} = zeros(pages,pages);

    for t=1:num_label_ERS(i)
        for j = 1:pages-1
            Vj = band_feature_ERS{i}{t}(:,:,j);
            Vj(isnan(Vj)) = [];
            for k = j+1:pages
                Vk =band_feature_ERS{i}{t}(:,:,k);
                Vk(isnan(Vk)) = [];
                sim_KL_ERS{i}(j,k,t) = KL(Vj(:)'+eps,Vk(:)'+eps);
                sim_L2_ERS{i}(j,k,t) = norm((Vj-Vk),2)/pages;
                sim_KL_ERS{i}(k,j,t) = sim_KL_ERS{i}(j,k,t);
                sim_L2_ERS{i}(k,j,t) = sim_L2_ERS{i}(j,k,t);
            end
        end
    end
    
    for j = 1:pages-1
        for k = j+1:pages
            temp2 = sim_L2_ERS{i}(j,k,:);
            sim_L21_ERS{i}(j,k) = norm(temp2(:),1);
            sim_L21_ERS{i}(k,j) =  sim_L21_ERS{i}(j,k);
        end
    end
    
    for j = 1:pages-1
        for k = j+1:pages
            temp2 = sim_KL_ERS{i}(j,k,:);
            sim_KL1_ERS{i}(j,k) = norm(temp2(:),1);
            sim_KL1_ERS{i}(k,j) =  sim_KL1_ERS{i}(j,k);
        end
    end
    
%     figure;
%     imagesc(sim_L21_ERS{i});
%     figure;
%     imagesc(sim_KL1_ERS{i}); 
end

save('Salinas_sim_L21_ERS','sim_L21_ERS');
save('Salinas_sim_KL1_ERS','sim_KL1_ERS');