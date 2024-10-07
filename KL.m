function D = KL(p, q)
    %KL Computes the Kullback-Leibler divergence between two distributions
    %   D = KL(p, q) computes the KL divergence of q from p.
    %
    %   Inputs:
    %       p - Probability distribution vector (must sum to 1)
    %       q - Probability distribution vector (must sum to 1)
    %
    %   Output:
    %       D - KL divergence value
    
        % Ensure p and q are column vectors
        p = p(:);
        q = q(:);
        
        % Add a small constant to avoid log(0) and division by zero
        epsilon = 1e-10;
        p = p + epsilon;
        q = q + epsilon;
        
        % Normalize the distributions to ensure they sum to 1
        p = p / sum(p);
        q = q / sum(q);
        
        % Compute KL divergence
        D = sum(p .* log(p ./ q));
    end