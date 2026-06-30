function [mua, iter] = DSI2mua_iterRecov(rhos, II, musp, opts)
% DSI2mua_iterRecov Recover absorption coefficient using dual-slope intensity iterative recovery (ZB boundary).
%
% [mua, iter] = DSI2mua_iterRecov(rhos, II, musp, opts)
%
% Written by Giles Blaney (Summer 2020; Ph.D. awarded May 2022)
%
% Inputs:
%   rhos - Source-detector distances in format [S1, L1, S2, L2] [mm]
%   II   - DC data in format [S1, L1, S2, L2] [mixed]
%   musp - Assumed reduced scattering coefficient [1/mm]
%   opts - (Optional) Structure containing recovery options [-]
%
% Outputs:
%   mua  - Recovered absorption coefficient [1/mm]
%   iter - Structure containing iteration information [-]

    %% Parse Input
    if nargin<=3
        mua0Bool=true;
        mueff_tol=1e-4; %1/mm
        n_max=10;
    else
        mua0Bool=false;
        mua0=opts.mua0;
        mueff_tol=opts.mueff_tol;
        n_max=opts.n_max;
    end
    
    if length(rhos)<4
        rhos=[rhos, rhos];
    end
    
    rho1=rhos(1:2);
    rho2=rhos(3:4);
    I1=II(1:2);
    I2=II(3:4);
    
    %% Find mua0
    if mua0Bool
        SSI10=(log(rho1(2)^2*I1(2))-log(rho1(1)^2*I1(1)))/...
            (rho1(2)-rho1(1));
        SSI20=(log(rho2(2)^2*I2(2))-log(rho2(1)^2*I2(1)))/...
            (rho2(2)-rho2(1));

        DSI0=(SSI10+SSI20)/2;

        mua0=DSI0^2/(3*musp);
    end
    
    %% Iterative mua Recovery
    z0=1/musp;
    x1=sqrt(z0^2+rho1.^2);
    x2=sqrt(z0^2+rho2.^2);

    n=1;
    stopBool=false;
    mueff0=sqrt(3*mua0*(musp+mua0));
    mueff=mueff0;
    while ~stopBool
        y1=log((I1.*x1.^2)./(mueff(n)+1./x1));
        y2=log((I2.*x2.^2)./(mueff(n)+1./x2));

        SSI1=diff(y1)/diff(x1);
        SSI2=diff(y2)/diff(x2);

        DSI=(SSI1+SSI2)/2;

        n=n+1;

        mueff(n)=-DSI;

        if n>=n_max
            warning('Max iterations reached.');
            stopBool=true;
        elseif abs(mueff(end)-mueff(end-1))<mueff_tol
            stopBool=true;
        end
        
        if n>=5 %Check for oscillation
            if sum(diff(mueff((end-2):end)))<=1e-10
                [~, ind]=min(abs(mueff((end-1):end)-mueff0));
                
                stopBool=true;
                
                mueff(end+1)=mueff(n+ind-2);
            end
        end
    end
    
    %% Package Output
    mua=sqrt(musp^2/4+mueff(end)^2/3)-musp/2;
    
    iter.n=n;
    iter.mueff_all=mueff;
    
end