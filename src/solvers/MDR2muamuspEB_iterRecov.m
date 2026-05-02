function [mua, musp, iter] = MDR2muamuspEB_iterRecov(rhos, R, opts)
% MDR2muamuspEB_iterRecov Recover absorption and reduced scattering using multi-distance reflectance iterative recovery (EB boundary).
%
% [mua, musp, iter] = MDR2muamuspEB_iterRecov(rhos, R, opts)
%
% Written by Giles Blaney, Ph.D. Fall 2020
%
% Inputs:
%   rhos - Vector of source-detector distances [mm]
%   R    - Vector of complex reflectance data [1/mm^2]
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
        nin=1.4;
        v=2.99792458e11/nin; %mm/sec
        nout=1;
    else
        mu0Bool=true;
%         mua0=opts.mua0;
%         musp0=opts.musp0;
        mueff_tol=opts.mueff_tol;
        n_max=opts.n_max;
        omega=opts.omega;
        v=opts.v;
        nin=opts.nin;
        nout=opts.nout;
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
    mua=mua0; %1/mm
    musp=musp0; %1/mm
    mueff0=sqrt(3*(mua0-1i*omega/v)*(musp0+mua0));
    A=n2A(nin, nout);
    mueff=mueff0;
    while ~stopBool
        z0=1/musp(n);
        zb=-2*A/(3*(musp(n)+mua(n)));
        z0p=-z0+2*zb;
        r1=sqrt(rhos.^2+z0.^2);
        r2=sqrt(rhos.^2+z0p.^2);

        C1=z0.*(1./r1+mueff(n))./(r1.^2);
        C2=-z0p.*(1./r2+mueff(n))./(r2.^2);

        x=r1;
        y=log(4*pi*R./...
            (C1+C2.*exp(mueff(n).*(r1-r2))));

        X=[x, ones(length(x), 1)];
        Y=pinv(X)*y;
        SR=Y(1);

        n=n+1;

        mueff(n)=-SR;

        a=real(mueff(n));
        b=imag(mueff(n));
        mutp=-2*a*b*v/(3*omega);
        mua(n)=(a^2-b^2)/(3*mutp);
        musp(n)=(mutp-mua(n));

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