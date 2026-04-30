function [mua, iter]=DSI2muaEB_iterRecov(rhos, II, musp, opts)
% [mua, iter]=DSI2muaEB_iterRecov(rhos, II, musp, opts)
% Giles Blaney Fall 2020
% Assumes extrapolated-boundary condition
% Expects DS set made of 4 SD measurements
%   Inputs:
%       dis      - 1 x 4 array of source-detector distances in the
%                  following format: [S1, L1, S2, L2] (mm)            
%       II       - DC data in a 1 x 4 array containing data in the
%                  following format: [S1, L1, S2, L2]
%       musp     - Assumed musp (1/mm)
%       opts     - (Optional) Structure containing options as feilds:
%                  - mua0: mua to use at for first iteration (1/mm)
%                  - mueff_tol: Change in mueff stoping criteria (1/mm)
%                    abs(mueff(n)-mueff(n-1))<mueff_tol => Stop
%                  - n_max: Maximum iteration number stoping criteria
%                    n>=n_max => Stop
%                  - nin: Index of refraction inside medium
%                  - nout: Index of refraction outside medium
%   Output:
%       mua      - mua found (1/mm)
%       iter     - Structure containing information about iterations as
%                  feilds:
%                  - n: Number of iterations
%                  - mueff_all: History of mueff during iterations

    %% Parse Input
    if nargin<=3
        mua0Bool=true;
        mueff_tol=1e-4; %1/mm
        n_max=10;
        nin=1.4;
        nout=1;
    else
        if isfield(opts, 'mua0')
            mua0=opts.mua0;
            mua0Bool=false;
        else
            mua0Bool=true;
        end
        mueff_tol=opts.mueff_tol;
        n_max=opts.n_max;
        nin=opts.nin;
        nout=opts.nout;
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
    r1_1=sqrt(z0^2+rho1.^2);
    r1_2=sqrt(z0^2+rho2.^2);

    n=1;
    stopBool=false;
    mua=mua0; %1/mm
    mueff0=sqrt(3*mua0*(musp+mua0));
    A=n2A(nin, nout);
    mueff=mueff0;
    while ~stopBool
        zb=-2*A/(3*(musp+mua(n)));
        z0p=-z0+2*zb;
        r2_1=sqrt(rho1.^2+z0p.^2);
        r2_2=sqrt(rho2.^2+z0p.^2);
        
        C1_1=z0.*(1./r1_1+mueff(n))./(r1_1.^2);
        C1_2=z0.*(1./r1_2+mueff(n))./(r1_2.^2);
        C2_1=-z0p.*(1./r2_1+mueff(n))./(r2_1.^2);
        C2_2=-z0p.*(1./r2_2+mueff(n))./(r2_2.^2);
        
        y_1=log(4*pi*I1./...
            (C1_1+C2_1.*exp(mueff(n).*(r1_1-r2_1))));
        y_2=log(4*pi*I2./...
            (C1_2+C2_2.*exp(mueff(n).*(r1_2-r2_2))));

        SSI1=diff(y_1)/diff(r1_1);
        SSI2=diff(y_2)/diff(r1_2);

        DSI=(SSI1+SSI2)/2;

        n=n+1;

        mueff(n)=-DSI;
        mua(n)=sqrt(musp^2/4+mueff(n)^2/3)-musp/2;
        
        if n>=n_max
%             warning('Max iterations reached.');
            stopBool=true;
        elseif abs(mueff(end)-mueff(end-1))<mueff_tol
            stopBool=true;
        end
        
        if n>=5 %Check for oscillation
            if sum(diff(mueff((end-2):end)))<=1e-10
                [~, ind]=min(abs(mueff((end-1):end)-mueff0));
                
                stopBool=true;
                
                mueff(end+1)=mueff(n+ind-2);
                mua(end+1)=mua(n+ind-2);
            end
        end
    end
    
    %% Package Output
    iter.n=n;
    iter.mueff_all=mueff;
    iter.mua=mua;
    
    mua=mua(end);
    
end