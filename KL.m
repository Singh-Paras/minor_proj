function D = vectorizedKL(P, Q)
    %VECTORIZEDKL Computes the KL divergence between multiple pairs of distributions
    %   D = VECTORIZEDKL(P, Q) computes the KL divergence for each pair of columns in P and Q
    %
    %   Inputs:
    %       P - Matrix where each column is a probability distribution
    %       Q - Matrix where each column is a probability distribution
    %
    %   Output:
    %       D - Vector of KL divergence values for each pair
    
        % Add epsilon to avoid log(0) and division by zero
        epsilon = 1e-10;
        P = P + epsilon;
        Q = Q + epsilon;
        
        % Normalize the distributions
        P = P ./ sum(P, 1);
        Q = Q ./ sum(Q, 1);
        
        % Compute KL divergence
        D = sum(P .* log(P ./ Q), 1);
    end