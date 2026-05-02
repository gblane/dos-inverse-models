function [mua, musp, iter]=MDR2muamusp_iterRecov(rhos, R, opts)
% [mua, musp, iter]=MDR2muamusp_iterRecov(rhos, R, opts)
% Giles Blaney Summer 2020
% Assumes zero-boundary condition (Equ. 12.A.7, Bigio & Fantini)
% For calibrated MD data with n (>=2) distances
%   Inputs:
%       rhos     - n x 1 array of source-detector distances (mm)            
%       R        - n x 1 array of complex relectance data
%       opts     - (Optional) Structure containing options as feilds:
%                  - mua0: mua to use at for first iteration (1/mm)
%                  - musp0: mua to use at for first iteration (1/mm)
%                  - mueff_tol: Change in mueff stoping criteria (1/mm)
%                    abs(mueff(n)-mueff(n-1))<mueff_tol => Stop
%                  - n_max: Maximum iteration number stoping criteria
%                    n>=n_max => Stop
%                  - omega: Modulation frequency (rad/sec) 
%                  - v: Speed of light in tissue (mm/sec)
%   Output:
%       mua      - mua found (1/mm)
%       mua      - musp found (1/mm)
%       iter     - Structure containing information about iterations as
%                  feilds:
%                  - n: Number of iterations
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
    
    if size(R, 1)==1
        R=R.';
    end
    if size(rhos, 1)==1
        rhos=rhos';
    end
    
    X=[rhos, ones(length(rhos), 1)];
    
    %% Find mua0
    if mu0Bool
        I=abs(R)/abs(R(1));
        lnr2I=log(rhos.^2.*I);
        P=angle(R)-angle(R(1));
        
        Y=X\lnr2I;
        SI0=Y(1);
        
        Y=X\P;
        SP0=Y(1);

        mua0=(omega/(2*v))*(SP0./SI0-SI0./SP0);
        musp0=(SI0.^2-SP0.^2)./...
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
        x=sqrt(z0^2+rhos.^2);
        X=[x, ones(length(x), 1)];
        
        y=log((R.*x.^2)./(z0*(mueff(n)+1./x)));
        
        Y=X\y;
        SR=Y(1);
        
        n=n+1;
        
        mueff(n)=-SR;
        
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