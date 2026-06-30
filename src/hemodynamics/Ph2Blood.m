function [dO, dD, dT, dmua, dPh, Ph0] = Ph2Blood(Ph, lambda, mua0, musp0, dis, bl, fmod, nin)
% Ph2Blood Calculate hemoglobin changes from phase measurements.
%
% [dO, dD, dT, dmua, dPh, Ph0] = Ph2Blood(Ph, lambda, mua0, musp0, dis, bl, fmod, nin)
%
% Written by Giles Blaney (Spring 2019; Ph.D. awarded May 2022)
%
% Inputs:
%   Ph     - Phase data [rad]
%   lambda - Optical wavelengths [nm]
%   mua0   - Baseline absorption coefficients [1/cm]
%   musp0  - Baseline reduced scattering coefficients [1/cm]
%   dis    - Source-detector distance [cm]
%   bl     - (Optional, default=1:end) Array of baseline indices [-]
%   fmod   - (Optional, default=1.40625e8 Hz) Modulation frequency [Hz]
%   nin    - (Optional, default=1.4) Internal index of refraction [-]
%
% Outputs:
%   dO     - Change in [HbO2] concentration [muM]
%   dD     - Change in [Hb] concentration [muM]
%   dT     - Change in [HbT] concentration [muM]
%   dmua   - Change in absorption coefficient [1/cm]
%   dPh    - Change in phase [rad]
%   Ph0    - Baseline phase [rad]
%
% Shared-repo dependencies:
%   circ_mean and wrapTo are provided by ../my-matlab.

    %% Setup
    if size(lambda, 1)~=1
        lambda=lambda';
    end

    if length(dis)~=length(lambda)
        dis=ones(size(lambda))*dis;
    elseif size(dis, 1)~=1
        dis=dis';
    end

    if size(Ph, 2)~=length(lambda)
        Ph=Ph';
    end
    
    if nargin<=5
        bl=1:size(Ph, 1);
    end
    if nargin<=6
        fmod=1.40625e8; %Hz
        nin=1.4;
    end
    L=complexPathLen(mean(dis), mua0, musp0, fmod, nin);
    Lph=imag(L);
    
    Ph=wrapTo((Ph),2*pi);
    Ph0=circ_mean(Ph(bl, :));
    Ph0=wrapTo(Ph0,2*pi);
    dPh=Ph-Ph0;
    
    %% mua
    dmua=-dPh./Lph;
    
    %% Blood
    E=makeE('OD', lambda)*1e4; % 1/(mM cm)
    X=linsolve(E, dmua');
    dO=X(1,:)'*1000; %uM
    dD=X(2,:)'*1000; %uM
    dT=dO+dD; %uM

end
