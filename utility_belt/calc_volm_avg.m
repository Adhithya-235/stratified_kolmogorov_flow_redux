function [datav] = calc_volm_avg(data,x,Lx,z,Lz)

    %   Performs volume averaging of provided data which is assumed to be a 3d array. 
    %   The trapezoidal rule is used, which is exact for periodic data provided the endpoints 
    %   are included. 
    
    %% CALCULATE VOLUME AVERAGE
    
    datav = squeeze(trapz(z, trapz(x, data, 2), 1)/(Lx*Lz));              
    
    end