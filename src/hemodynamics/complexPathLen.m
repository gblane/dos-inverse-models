function [L] = complexPathLen(dis, mua0, musp0, fmod, nin)
% complexPathLen Calculate the complex pathlength for a semi-infinite medium.
%
% [L] = complexPathLen(dis, mua0, musp0, fmod, nin)
%
% Written by Giles Blaney, Ph.D. Spring 2019
%
% Inputs:
%   dis   - Source-detector distance [cm]
%   mua0  - Baseline absorption coefficients [1/cm]
%   musp0 - Baseline reduced scattering coefficients [1/cm]
%   fmod  - (Optional, default=1.40625e8 Hz) Modulation frequency [Hz]
%   nin   - (Optional, default=1.4) Internal index of refraction [-]
%
% Outputs:
%   L     - Complex pathlength [cm]

    warning('Legacy function DO NOT USE, use complexTotPathLen() instead');

    %% Setup
    if nargin<=3
        fmod=1.40625e8; %Hz
        nin=1.4;
    end
    nout=1;
    c=2.99792458e10; %cm/sec
    A=n2A(nin, nout);
    omega=2*pi*fmod; %rad/sec
    nu=c/nin;
    
    %% Calc <L>
    L=ones(size(mua0))*NaN;
    for lamInd=1:length(mua0)
        mua=mua0(lamInd);
        musp=musp0(lamInd);
        
        x0=1/musp; %cm
        D=1/(3*musp); %cm
        xb=-2*A*D; %cm

        mueff=sqrt(mua/D-1i*omega/(nu*D)); %1/cm

        rs=[x0, 0]; %cm
        rsp=[-x0+2*xb, 0]; %cm

        ri=[0, dis]; %cm

        r1=norm(ri-rs);
        r2=norm(ri-rsp);

        R=(x0*(1/r1+mueff)*(exp(-mueff*r1)/(r1^2))+...
            (x0-2*xb)*(1/r2+mueff)*(exp(-mueff*r2)/(r2^2)))/(4*pi); %1/cm^2

        L(lamInd)=((x0/r1)*exp(-mueff*r1)+...
            ((x0-2*xb)/r2)*exp(-mueff*r2))/(8*pi*D*R); %cm
    end
end