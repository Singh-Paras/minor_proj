clc; clear; close all;

% Load data
load Salinas_corrected.mat;
load Salinas_pca_3.mat;

% Preprocess data
cube_data = double(salinas_corrected);
[rows, cols, pages] = size(cube_data);
cube_data_flat = reshape(cube_data, rows * cols, pages)';    % Reshape 3D to 2D
cube_data_nm_flat = mapminmax(cube_data_flat, 0, 255);      % Normalize data to [0, 255]
cube_data_nm = reshape(cube_data_nm_flat', rows, cols, pages);

% Parameters
T = 3.7; % Salinas  
K = [1, 2, 3];
Nk = floor(K .* (rows * cols) / (T * T));
sp_label_ERS = zeros(rows, cols, length(K));
num_label_ERS = zeros(1, length(K));

% Preallocate cell arrays
band_feature_ERS = cell(1, length(K));
sim_KL_ERS = cell(1, length(K));
sim_L2_ERS = cell(1, length(K));
sim_L21_ERS = cell(1, length(K));
sim_KL1_ERS = cell(1, length(K));

%% Superpixel Segmentation ERS
for i = 1:length(K)
    sp_label_ERS(:, :, i) = mymex_ers(double(Salinas_pca_3), Nk(i), 0.05, 5);  % ERS superpixel segmentation
    num_label_ERS(i) = max(sp_label_ERS(:, :, i), [], 'all') + 1;
    
    %% Compute Distances Between Superpixels
    band_feature_ERS{i} = SpRectangle(sp_label_ERS(:, :, i), num_label_ERS(i), cube_data_nm);
    sim_KL_ERS{i} = zeros(pages, pages, num_label_ERS(i));
    sim_L2_ERS{i} = zeros(pages, pages, num_label_ERS(i));
    sim_L21_ERS{i} = zeros(pages, pages);
    sim_KL1_ERS{i} = zeros(pages, pages);
    
    % Preallocate temporary storage for vectorized KL computation
    temp_P = [];
    temp_Q = [];
    
    for t = 1:num_label_ERS(i)
        for j = 1:pages-1
            Vj = band_feature_ERS{i}{t}(:, :, j);
            Vj(isnan(Vj)) = [];
            for k = j+1:pages
                Vk = band_feature_ERS{i}{t}(:, :, k);
                Vk(isnan(Vk)) = [];
                
                if isempty(Vj) || isempty(Vk)
                    continue;
                end
                
                % Accumulate vectors
                temp_P = [temp_P; Vj(:)' + eps];
                temp_Q = [temp_Q; Vk(:)' + eps];
            end
        end
    end
    
    % Compute KL divergence in bulk
    if ~isempty(temp_P) && ~isempty(temp_Q)
        D_KL = vectorizedKL(temp_P', temp_Q');
    else
        D_KL = [];
    end
    
    % Assign the computed KL divergence back to the similarity matrix
    idx = 1;
    for t = 1:num_label_ERS(i)
        for j = 1:pages-1
            for k = j+1:pages
                if idx > length(D_KL)
                    break;
                end
                sim_KL_ERS{i}(j, k, t) = D_KL(idx);
                sim_KL_ERS{i}(k, j, t) = D_KL(idx);
                sim_L2_ERS{i}(j, k, t) = norm(band_feature_ERS{i}{t}(:, :, j) - band_feature_ERS{i}{t}(:, :, k), 2) / pages;
                sim_L2_ERS{i}(k, j, t) = sim_L2_ERS{i}(j, k, t);
                idx = idx + 1;
            end
        end
    end
    
    % Compute aggregated similarities
    for j = 1:pages-1
        for k = j+1:pages
            temp2 = sim_L2_ERS{i}(j, k, :);
            sim_L21_ERS{i}(j, k) = norm(temp2(:), 1);
            sim_L21_ERS{i}(k, j) = sim_L21_ERS{i}(j, k);
        end
    end
    
    for j = 1:pages-1
        for k = j+1:pages
            temp2 = sim_KL_ERS{i}(j, k, :);
            sim_KL1_ERS{i}(j, k) = norm(temp2(:), 1);
            sim_KL1_ERS{i}(k, j) = sim_KL1_ERS{i}(j, k);
        end
    end
    
    % Clear temporary variables to save memory
    clear temp_P temp_Q D_KL
end

% Save the results
save('Salinas_sim_L21_ERS.mat', 'sim_L21_ERS');
save('Salinas_sim_KL1_ERS.mat', 'sim_KL1_ERS');