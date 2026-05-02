function [mua, musp, iter] = DSR2muamuspEB_iterRecov(rhos, RR, opts)
% DSR2muamuspEB_iterRecov Recover absorption and reduced scattering using dual-slope reflectance iterative recovery (EB boundary).
%
% [mua, musp, iter] = DSR2muamuspEB_iterRecov(rhos, RR, opts)
%
% Written by Giles Blaney, Ph.D. Fall 2020
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
%                  feilds:
%                  - n: Number of iterations
%                  - mueff_all: History of mueff during iterations
%                  - mua: History of mua during iterations
%                  - musp: History of musp during iterations

    %% Parse Input
    if nargin<=2
        mu0Bool=true;
        mueff_tol=1e-4; %1/mm
        n_max=10;
        omega=2*pi*140.625e6; %rad/sec
        nin=1.4;
        v=2.99792458e11/nin; %mm/sec
        nout=1;
        enforcePosMUA=false;
    else
        if isfield(opts, 'mua0')
            mua0=opts.mua0;
            musp0=opts.musp0;
            mu0Bool=false;
        else
            mu0Bool=true;
        end
        if isfield(opts, 'mueff_tol')
            mueff_tol=opts.mueff_tol;
        else
            mueff_tol=1e-4; %1/mm
        end
        if isfield(opts, 'n_max')
            n_max=opts.n_max;
        else
            n_max=10;
        end
        if isfield(opts, 'omega')
            omega=opts.omega;
        else
            omega=2*pi*140.625e6; %rad/sec
        end
        if isfield(opts, 'nin')
            nin=opts.nin;
        else
            nin=1.4;
        end
        if isfield(opts, 'nout')
            nout=opts.nout;
        else
            nout=1;
        end
        if isfield(opts, 'v')
            v=opts.v;
        else
            v=2.99792458e11/nin; %mm/sec
        end
        if isfield(opts, 'enforcePosMUA')
            enforcePosMUA=opts.enforcePosMUA;
        else
            enforcePosMUA=false;
        end
    end
    
    if length(rhos)<4
        rhos=[rhos, rhos];
    end
    
    rho1=rhos(1:2);
    rho2=rhos(3:4);
    R1=RR(1:2);
    R2=RR(3:4);
    
    if angle(R1(1))>angle(R1(2))
        R1(2)=abs(R1(2))*exp(1i*(angle(R1(2))+2*pi));        
    end
        
    if angle(R2(1))>angle(R2(2))
        R2(2)=abs(R2(2))*exp(1i*(angle(R2(2))+2*pi));        
    end
    
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
        
        SSP10=wrapToPi(P1(2)-P1(1))/...
            (rho1(2)-rho1(1));
        SSP20=wrapToPi(P2(2)-P2(1))/...
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
    mua=mua0; %1/mm
    musp=musp0; %1/mm
    mueff0=sqrt(3*(mua0-1i*omega/v)*(musp0+mua0));
    A=n2A(nin, nout);
    mueff=mueff0;
    while ~stopBool
        z0=1/musp(n);
        zb=-2*A/(3*(musp(n)+mua(n)));
        z0p=-z0+2*zb;
        r1_1=sqrt(rho1.^2+z0.^2);
        r1_2=sqrt(rho2.^2+z0.^2);
        r2_1=sqrt(rho1.^2+z0p.^2);
        r2_2=sqrt(rho2.^2+z0p.^2);
        
        C1_1=z0.*(1./r1_1+mueff(n))./(r1_1.^2);
        C1_2=z0.*(1./r1_2+mueff(n))./(r1_2.^2);
        C2_1=-z0p.*(1./r2_1+mueff(n))./(r2_1.^2);
        C2_2=-z0p.*(1./r2_2+mueff(n))./(r2_2.^2);
        
        y_1=log(4*pi*R1./...
            (C1_1+C2_1.*exp(mueff(n).*(r1_1-r2_1))));
        y_2=log(4*pi*R2./...
            (C1_2+C2_2.*exp(mueff(n).*(r1_2-r2_2))));

        SSR1=diff(y_1)/diff(r1_1);
        SSR2=diff(y_2)/diff(r1_2);
        
        DSR=(SSR1+SSR2)/2;

        n=n+1;

        mueff(n)=-DSR;
        
        a=real(mueff(n));
        b=imag(mueff(n));
        mutp=-2*a*b*v/(3*omega);
        mua(n)=(a^2-b^2)/(3*mutp);
        if enforcePosMUA
            mua(n)=max([0, mua(n)]);
        end
        musp(n)=(mutp-mua(n));
        
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