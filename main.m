clc; clear; close all; warning('off');

%% 加载高光谱数据

load Salinas_corrected.mat;
cube_data = double(salinas_corrected);
[rows, cols, pages] = size(cube_data);
cube_data_flat = reshape(cube_data, rows * cols, pages)';	% 三维变二维
cube_data_nm_flat = mapminmax(cube_data_flat, 0, 255) ;  % 数据归一化0~255
cube_data_nm = reshape(cube_data_nm_flat', rows, cols, pages);

K = [1,2,3];
Nk = zeros(1,length(K));

load Salinas_sim_L21_ERS.mat;
load Salinas_sim_KL1_ERS.mat;
load Salinas_sim_L21_SNIC.mat;
load Salinas_sim_KL1_SNIC.mat;
load Salinas_sim_L21_SLIC.mat;
load Salinas_sim_KL1_SLIC.mat;

%% 分组
group_ers = cell(1,length(Nk));
band_idx_sep_ers = cell(1,length(Nk));
group_snic = cell(1,length(Nk));
band_idx_sep_snic = cell(1,length(Nk));
group_slic = cell(1,length(Nk));
band_idx_sep_slic = cell(1,length(Nk));


for i = 1:length(Nk)
    [locs_ERS{i}, ~] = sub_part(sim_KL1_ERS{i},pages); 
    group_ers{i} = [1,locs_ERS{i},pages];    

    [locs_SNIC{i}, ~] = sub_part(sim_KL1_SNIC{i},pages); 
    group_snic{i} = [1,locs_SNIC{i},pages];   

    [locs_SLIC{i}, ~] = sub_part(sim_KL1_SLIC{i},pages); 
    group_slic{i} = [1,locs_SLIC{i},pages];
end

%% 聚类

ent = Entrop(cube_data_nm_flat');  % 信息熵
for i = 1:length(Nk)  
% ERS
    for j = 1:length(group_ers{i})-1
        if j == 1
            sim_sep_ERS = sim_L21_ERS{i}(group_ers{i}(j):group_ers{i}(j+1),group_ers{i}(j):group_ers{i}(j+1));
            band_idx_sep_ers{i} = [band_idx_sep_ers{i}, sufdpc(sim_sep_ERS,ent(group_ers{i}(j):group_ers{i}(j+1)),pages)+group_ers{i}(j)-1];
        else
            sim_sep_ERS = sim_L21_ERS{i}(group_ers{i}(j)+1:group_ers{i}(j+1),group_ers{i}(j)+1:group_ers{i}(j+1));
            band_idx_sep_ers{i} = [band_idx_sep_ers{i}, sufdpc(sim_sep_ERS,ent(group_ers{i}(j)+1:group_ers{i}(j+1)),pages)+group_ers{i}(j)];
        end
    end
% SNIC
    for j = 1:length(group_snic{i})-1
        if j == 1
            sim_sep_SNIC = sim_L21_SNIC{i}(group_snic{i}(j):group_snic{i}(j+1),group_snic{i}(j):group_snic{i}(j+1));
            band_idx_sep_snic{i} = [band_idx_sep_snic{i}, sufdpc(sim_sep_SNIC,ent(group_snic{i}(j):group_snic{i}(j+1)),pages)+group_snic{i}(j)-1];
        else
            sim_sep_SNIC = sim_L21_SNIC{i}(group_snic{i}(j)+1:group_snic{i}(j+1),group_snic{i}(j)+1:group_snic{i}(j+1));
            band_idx_sep_snic{i} = [band_idx_sep_snic{i}, sufdpc(sim_sep_SNIC,ent(group_snic{i}(j)+1:group_snic{i}(j+1)),pages)+group_snic{i}(j)];
        end
    end
% SLIC
    for j = 1:length(group_slic{i})-1
        if j == 1
            sim_sep_SLIC = sim_L21_SLIC{i}(group_slic{i}(j):group_slic{i}(j+1),group_slic{i}(j):group_slic{i}(j+1));
            band_idx_sep_slic{i} = [band_idx_sep_slic{i}, sufdpc(sim_sep_SLIC,ent(group_slic{i}(j):group_slic{i}(j+1)),pages)+group_slic{i}(j)-1];
        else
            sim_sep_SLIC = sim_L21_SLIC{i}(group_slic{i}(j)+1:group_slic{i}(j+1),group_slic{i}(j)+1:group_slic{i}(j+1));
            band_idx_sep_slic{i} = [band_idx_sep_slic{i}, sufdpc(sim_sep_SLIC,ent(group_slic{i}(j)+1:group_slic{i}(j+1)),pages)+group_slic{i}(j)];
        end
    end
end

%% 投票

band_num = 10; 

band_idx_sepp = cell2mat(band_idx_sep_ers);
band_idx_sepp = [band_idx_sepp, cell2mat(band_idx_sep_snic)];
band_idx_sepp = [band_idx_sepp, cell2mat(band_idx_sep_slic)];
temp3 = tabulate(band_idx_sepp(:));
temp3 = sortrows(temp3,size(temp3,2),'descend');
if(~isempty(temp3))
    score = ent(temp3(find(temp3(:,2)~=0),1)).*temp3(find(temp3(:,2)~=0),2);
    [~, sortIndex2] = sort(score(:), 'descend');
    maxIndex2 = sortIndex2(1:band_num);
    msgcf_bs = temp3(maxIndex2)';
end
