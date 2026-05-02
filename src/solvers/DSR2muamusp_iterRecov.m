function [mua, musp, iter] = DSR2muamusp_iterRecov(rhos, RR, opts)
% DSR2muamusp_iterRecov Recover absorption and reduced scattering using dual-slope reflectance iterative recovery (ZB boundary).
%
% [mua, musp, iter] = DSR2muamusp_iterRecov(rhos, RR, opts)
%
% Written by Giles Blaney, Ph.D. Summer 2020
%
% Inputs:
%   rhos - Source-detector distances in format [S1, L1, S2, L2] [mm]
%   RR   - Complex reflectance data in format [S1, L1, S2, L2] [1/mm^2]
%   opts - (Optional) Structure containing recovery options [-]
%
% Outputs:
%   mua  - Recovered absorption coefficient [1/mm]
%   musp - Recovered reduced scattering coefficient [1/mm]
%   iter - Structure containing iteration information [-]
%                  - mueff_all: History of mueff during iterations
%                  - mua_all: History of mua during iterations
%                  - musp_all: History of musp during iterations

    %% Parse Input
    if nargin<=2
        mu0Bool=true;
        mueff_tol=1e-4; %1/mm
        n_max=10;
        omega=2*pi*140.625e6; %rad/sec
        v=2.99792458e11/1.4; %mm/sec
    else
        mu0Bool=false;
        mua0=opts.mua0;
        musp0=opts.musp0;
        mueff_tol=opts.mueff_tol;
        n_max=opts.n_max;
        omega=opts.omega;
        v=opts.v;
    end
    
    if length(rhos)<4
        rhos=[rhos, rhos];
    end
    
    rho1=rhos(1:2);
    rho2=rhos(3:4);
    R1=RR(1:2);
    R2=RR(3:4);
    
    %% Find mua0
    if mu0Bool
        I1=abs(R1);
        I2=abs(R2);
        P1=angle(R1);
        P2=angle(R2);
        
        SSI10=(log(rho1(2)^2*I1(2))-log(rho1(1)^2*I1(1)))/...
            (rho1(2)-rho1(1));
        SSI20=(log(rho2(2)^2*I2(2))-log(rho2(1)^2*I2(1)))/...
            (rho2(2)-rho2(1));
        
        SSP10=(P1(2)-P1(1))/...
            (rho1(2)-rho1(1));
        SSP20=(P2(2)-P2(1))/...
            (rho2(2)-rho2(1));

        DSI0=(SSI10+SSI20)/2;
        DSP0=(SSP10+SSP20)/2;

        mua0=(omega/(2*v))*(DSP0./DSI0-DSI0./DSP0);
        musp0=(DSI0.^2-DSP0.^2)./...
            (3*mua0)-mua0;
    end
    
    %% Iterative mua Recovery
    n=1;
    stopBool=false;
    mua=mua0;
    musp=musp0;
    mueff0=sqrt(3*(mua0-1i*omega/v)*(musp0+mua0));
    mueff=mueff0;
    z0=1/musp;
    while ~stopBool
        x1=sqrt(z0^2+rho1.^2);
        x2=sqrt(z0^2+rho2.^2);
        
        y1=log((R1.*x1.^2)./(z0*(mueff(n)+1./x1)));
        y2=log((R2.*x2.^2)./(z0*(mueff(n)+1./x2)));

        SSR1=diff(y1)/diff(x1);
        SSR2=diff(y2)/diff(x2);

        DSR=(SSR1+SSR2)/2;

        n=n+1;

        mueff(n)=-DSR;
        
        a=real(mueff(n));
        b=imag(mueff(n));
        mutp=-2*a*b*v/(3*omega);
        mua(n)=(a^2-b^2)/(3*mutp);
        musp(n)=(mutp-mua(n));
        z0=1/musp(n);

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
                
                if musp(n+ind-2)>0
                    mua(end+1)=mua(n+ind-2);
                    musp(end+1)=musp(n+ind-2);
                else
                    ind=mod(ind, 2)+1;
                    if musp(n+ind-2)>0
                        mua(end+1)=mua(n+ind-2);
                        musp(end+1)=musp(n+ind-2);
                    else
                        mua(end+1)=NaN;
                        musp(end+1)=NaN;
                    end
                end
            end
        end
    end
    
    %% Package Output
    iter.n=n;
    iter.mueff_all=mueff;
    iter.mua=mua;
    iter.musp=musp;
    
    mua=mua(end);
    musp=musp(end);
    
end