clc; clear; close all; warning('off');

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
numK = length(K);
Nk = floor(K .* (rows * cols) / (T * T));

% Preallocate cell arrays
sp_label_ERS = zeros(rows, cols, numK);
num_label_ERS = zeros(1, numK);

band_feature_ERS = cell(1, numK);
sim_KL_ERS = cell(1, numK);
sim_L2_ERS = cell(1, numK);
sim_L21_ERS = cell(1, numK);
sim_KL1_ERS = cell(1, numK);

% Start Parallel Pool (if not already started)
pool = gcp('nocreate');
if isempty(pool)
    numWorkers = 4; % Adjust based on your system
    pool = parpool(numWorkers);
    fprintf('Started a parallel pool with %d workers.\n', numWorkers);
else
    fprintf('Parallel pool already running with %d workers.\n', pool.NumWorkers);
end

%% Superpixel Segmentation ERS using parfor
parfor i = 1:numK
    fprintf('Processing K = %d (%d of %d)\n', K(i), i, numK);
    
    % Superpixel segmentation
    sp_label = mymex_ers(double(Salinas_pca_3), Nk(i), 0.05, 5);  % ERS superpixel segmentation
    sp_label_ERS(:,:,i) = sp_label;
    num_labels = max(sp_label, [], 'all') + 1;
    num_label_ERS(i) = num_labels;
    
    % Compute Distances Between Superpixels
    band_feature = SpRectangle(sp_label, num_labels, cube_data_nm);
    band_feature_ERS{i} = band_feature;
    
    sim_KL = zeros(pages, pages, num_labels);
    sim_L2 = zeros(pages, pages, num_labels);
    
    for t = 1:num_labels
        for j = 1:pages-1
            Vj = band_feature{t}(:, :, j);
            Vj(isnan(Vj)) = [];
            for k = j+1:pages
                Vk = band_feature{t}(:, :, k);
                Vk(isnan(Vk)) = [];
                
                if isempty(Vj) || isempty(Vk)
                    continue;
                end
                
                % Compute KL divergence
                sim_KL(j, k, t) = KL(Vj(:)' + eps, Vk(:)' + eps);
                sim_L2(j, k, t) = norm(Vj - Vk, 2) / pages;
                
                % Assign symmetric values
                sim_KL(k, j, t) = sim_KL(j, k, t);
                sim_L2(k, j, t) = sim_L2(j, k, t);
            end
        end
    end
    
    sim_KL_ERS{i} = sim_KL;
    sim_L2_ERS{i} = sim_L2;
    
    % Compute aggregated similarities
    sim_L21 = zeros(pages, pages);
    sim_KL1 = zeros(pages, pages);
    
    for j = 1:pages-1
        for k = j+1:pages
            temp_L2 = sim_L2(j, k, :);
            sim_L21(j, k) = norm(temp_L2(:), 1);
            sim_L21(k, j) = sim_L21(j, k);
            
            temp_KL = sim_KL(j, k, :);
            sim_KL1(j, k) = norm(temp_KL(:), 1);
            sim_KL1(k, j) = sim_KL1(j, k);
        end
    end
    
    sim_L21_ERS{i} = sim_L21;
    sim_KL1_ERS{i} = sim_KL1;
    
    fprintf('Completed processing for K = %d (%d of %d)\n', K(i), i, numK);
end

% Save the results
save('Salinas_sim_L21_ERS.mat', 'sim_L21_ERS');
save('Salinas_sim_KL1_ERS.mat', 'sim_KL1_ERS');

% Close the parallel pool (optional)
delete(pool);