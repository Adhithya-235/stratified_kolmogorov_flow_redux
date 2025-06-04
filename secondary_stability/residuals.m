function T = residuals(eigvals, eigvecs, Afun)
%COMPUTE_EIGEN_RESIDUALS Compute residual norms for a set of eigenpairs
%
%   T = compute_eigen_residuals(eigvals, eigvecs, Afun)
%
%   Inputs:
%     eigvals : k x 1 vector of eigenvalues (complex or real)
%     eigvecs : n x k matrix whose columns are eigenvectors
%     Afun    : function handle such that Afun(v) = A*v
%
%   Output:
%     T : table containing eigenvalue index, real part, imaginary part,
%         and residual norm for each eigenpair

    %% Number of eigenpairs
    
    k = length(eigvals);

    %% Preallocate residual array
    
    residuals = zeros(k, 1);

    %% Compute residuals
    
    for j = 1:k
        lambda = eigvals(j);
        v = eigvecs(:, j);
        Av = Afun(v);
        residuals(j) = norm(Av - lambda * v)/abs(lambda);
    end

    %% Create output table
    
    T = table((1:k)', real(eigvals), imag(eigvals), residuals, ...
        'VariableNames', {'Index', 'RealPart', 'ImagPart', 'ResidualNorm'});

    % Display the table
    
    disp(T);

end
