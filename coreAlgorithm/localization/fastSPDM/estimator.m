function [maxQ,mx,my,dx,dy,sx2,sy2,Q]=estimator(data,meanbg, SPDM_EMCCD_Variance)
% function [maxQ,mx,my,dx,dy,sx,sy,Q]=estimator(data,meanbg)
%   The Gaussian estimator calculates the total intensity Q, the peak intensity maxQ,
%   the center-of-mass mx, my, the confidence (fit error) dy, dy, and the squared width
%   of the distribution sx2, sy2
%
%   For encreased performance the Matlab function sum() and dot products are used instead of for loops
% 
%   to calculate the sigma and error the original data is used.
%   
%   Output:
%       maxQ    -maximum of gauss function
%       mx,my   -estimated meanvalues
%       dx,dy   -error of meanvalues
%       sx,sy   -width of the distribution
%       Qs      -accumulated charge of roi
%
%created by: Manfred Kirchgessner < Manfred.Kirchgessner@web.de>, 
%            Frederik Grüll <Frederik.Gruell@kip.uni-heidelberg.de>%
%


    [Sizey Sizex] = size(data);

    data = data(:);

    %% Qges
    Q = sum(data);

    persistent xa;
    persistent ya;
   
    if (isempty(xa))
        [xa, ya] = meshgrid( 1:Sizex, 1:Sizey);

        xa=xa(:);
        ya=ya(:);
        
    end
    
    %% Accumulators
    qx = data.*xa;
    qx2_acc =sum((qx.*xa));
    qx_acc = sum(qx);
    
    qy = data.*ya;
    qy2_acc = sum(qy.*ya);
    qy_acc = sum(qy);
    
    q_gr_0 = (data>0);

    x = xa.*q_gr_0;
    y = ya.*q_gr_0;

    x2_acc = sum(xa.*x);
    y2_acc = sum(ya.*y);
    x_acc = sum(x);
    y_acc = sum(y);
    
    n_acc = sum(q_gr_0);
    
    %% mx my
    mx = qx_acc / Q;
    my = qy_acc / Q;
    %% sx2 sy2
    sx2 = max(qx2_acc / Q - mx * mx,1/12);
    sy2 = max(qy2_acc / Q - my * my,1/12);

    Q2 = Q*Q;
    
    %% dx dy
    %localization precision according to Thompson et al.
    %dx = SPDM_EMCCD_Variance*sqrt(  (1/12 + sx2)/Q + ...
    %            (x2_acc - mx *(2 * x_acc - n_acc *mx) ) * meanbg / Q2  );
    %dy = SPDM_EMCCD_Variance*sqrt(  (1/12 + sy2)/Q + ...
    %            (y2_acc - my *(2 * y_acc - n_acc *my) ) * meanbg / Q2  );
            
    %explicit formulation according to Thompson et al.
    %sigmax2 = (1/12 + sx2);
    %sigmay2 = (1/12 + sy2);
    %dx = SPDM_EMCCD_Variance*sqrt(sigmax2/Q * (1 + 8*pi*sigmax2*meanbg/Q));
    %dy = SPDM_EMCCD_Variance*sqrt(sigmay2/Q * (1 + 8*pi*sigmay2*meanbg/Q));
    
    %localization precision according to Mortensen et al. (approximation according to Stallinga et al.)
    taux = 2*pi*(sx2+1/12)*meanbg/Q;
    tauy = 2*pi*(sy2+1/12)*meanbg/Q;
    dx = SPDM_EMCCD_Variance*sqrt( (sx2+1/12)/Q * (1+4*taux+sqrt(2*taux/(1+4*taux))) );
    dy = SPDM_EMCCD_Variance*sqrt( (sy2+1/12)/Q * (1+4*tauy+sqrt(2*tauy/(1+4*tauy))) );
    %% Qmax
    maxQ = max(data);
    %maxQ = Q/(2*pi*sqrt(sx2*sy2));

end
