function [Xi, B, Psi] = get_ecs_fields(Reb, rescale)

    % Construct nonlinear basic state vector in grid space. 
    
    %% LOAD DATA FROM FILE
    
    datafile = sprintf('../import/ECS_real_field_for_Reb=%d.mat', Reb);
    nlbs     = load(datafile);
    
    %% GET ORIGINAL Nx, Nz
    
    Nx = length(nlbs.x);
    Nz = length(nlbs.z);
    
    %% PHYSICAL SPACE CONSTRUCTION
    
    Xi  = interpft(interpft(nlbs.omega, Nz*rescale, 1), Nx*rescale, 2);
    Psi = interpft(interpft(nlbs.psi, Nz*rescale, 1), Nx*rescale, 2);
    B   = interpft(interpft(nlbs.Buoyancy, Nz*rescale, 1), Nx*rescale, 2); 
    
end