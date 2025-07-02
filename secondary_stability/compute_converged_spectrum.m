function [eigvals, eigvecs, dom_mode, meta] = compute_converged_spectrum(...
    Reb, Pr, Fr, alpha, beta, numEig, residualtol, maxit, targetRelErr, initialRescale)

%% DEFAULT VALUES

if nargin < 10
    targetRelErr = 1e-4;
end
if nargin < 11
    initialRescale = 0.5;
end

%% INITIALIZE

relErr     = Inf;
vecErr     = Inf;
scale      = initialRescale;
prevEig    = NaN;
outputflag = 1;  % Verbose mode

%% MAIN LOOP

while (relErr > targetRelErr) || (vecErr > 1e-14)
    
    %% CHECK RESOLUTION CAP AND BREAK LOOP IF EXCEEDED
    
    if exist('meta','var') && any(meta.resolution >= 512)
        warning("Maximum allowed grid size (512 x 512) exceeded. Terminating refinement.");
        break;
    end

    fprintf('Solving ECS stability with Re_b = %.2f, alpha = %.2f, beta = %.2f, Scale = %.2f...\n', ...
        Reb, alpha, beta, scale);

    %% SOLVE EIGENVALUE PROBLEM

    [eigvals, eigvecs, dom_mode, meta] = stability_eigensolve(Reb, Pr, Fr, alpha, beta, ...
        numEig, residualtol, maxit, outputflag, scale);

    newEig = eigvals(1);  % dominant eigenvalue

    %% EIGENVALUE RELATIVE ERROR

    if ~isnan(prevEig)
        relErr = abs((newEig - prevEig) / prevEig);
        fprintf("Relative Error (dominant eigenvalue): %.3e\n", relErr);

        if relErr <= targetRelErr
            
            %% CHECK EIGENVECTOR ALIGNMENT
            
            Nx_prev = meta_prev.resolution(1); Nz_prev = meta_prev.resolution(2);
            Nx_curr = meta.resolution(1); Nz_curr = meta.resolution(2);
            vecErr = compare_eigenvectors_fftstack(eigvecs_prev, eigvecs, ...
                Nx_prev, Nz_prev, Nx_curr, Nz_curr);
            fprintf("Eigenvector mismatch (1 - cos θ): %.2e\n", vecErr);
            
            %% ITERATION RESULT LOGGING 
            
            if vecErr < 1e-12
                fprintf("✅ Eigenvector and eigenvalue convergence satisfied.\n");
            else
                warning("Eigenvectors are not parallel to machine precision — refining further.");
            end
        
        end
    
    end

    %% UPDATE FOR NEXT ITERATION

    prevEig = newEig;
    eigvecs_prev = eigvecs;
    meta_prev = meta;
    scale = scale * 2;

end

fprintf("Grid-independent spectrum found at scale: %.2f\n", scale/2);

end