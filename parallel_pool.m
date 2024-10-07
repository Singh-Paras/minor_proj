% Check if a parallel pool is already running
pool = gcp('nocreate'); % Get current pool without creating a new one

if isempty(pool)
    numWorkers = 4; % Adjust based on your system's CPU cores
    pool = parpool(numWorkers);
    fprintf('Started a parallel pool with %d workers.\n', numWorkers);
else
    fprintf('Parallel pool already running with %d workers.\n', pool.NumWorkers);
end