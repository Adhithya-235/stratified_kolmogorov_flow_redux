function OptScale = select_optimal_grid(Reb, Pr, Fr, alpha, beta,...
    residualtol, maxit, targetRelErr, initialRescale)

%% DEFAULT VALUES

if nargin < 8
    targetRelErr = 1e-4; 
end
if nargin < 9  
    initialRescale = 1/2; 
end

%% PARAMETERS

numEig     = 1;
outputflag = 0;
relErr     = Inf;
scale      = initialRescale;
prevEig    = NaN;

%% CONDUCT GRID INDEPENDENCE TESTS

while relErr > targetRelErr

    fprintf('Solving eigenvalue problem with Re_b = %f, Scale = %1.2f.\n',...
        Reb, scale)
    
    [eigvals, ~, ~, ~] = stability_eigensolve(Reb, Pr, Fr, alpha, beta,...
        numEig, residualtol, maxit, outputflag, scale);
    newEig = eigvals(1);  % dominant mode
    
    if ~isnan(prevEig)
        relErr = abs((newEig - prevEig) / prevEig);
        fprintf("Scale %.2f → Relative Error: %.3e\n", scale, relErr);
    end
    
    prevEig = newEig;
    scale = scale * 2;

end

fprintf("Grid-independent eigenvalue found: %.5e\n", newEig);

%% OUTPUT

OptScale = scale/2;

end