function [dO, dD, dT, dmua, dPh, Ph0] = Ph2Blood(Ph, lambda, mua0, musp0,...
    dis, bl, fmod, nin)
% Giles Blaney Spring 2019
% [dO, dD, dT, dmua, dI, I0] = Ph2Blood(Ph, lambda, dis, baseline)
% Inputs:
%   Ph     - Phase data time x wavelength array. (rad)
%   lambda - Optical wavelengths 1 x wavelength array. (nm)
%   mua0   - Baseline absorption in a 1 x wavelength array. (1/cm)
%   musp0  - Baseline scattering in a 1 x wavelength array. (1/cm)
%   dis    - Source detector distance. (cm)
%   bl     - (OPTIONAL, default=1:end) Array of baseline indexes.
%   fmod   - (OPTIONAL, default=1.40625e8 Hz) Modulation frequecy. (Hz)
%   nin    - (OPTIONAL, default=1.4) Internal index of refraction. (-)
% Outputs:
%   dO     - Change in [HbO2] concentration in a time x 1 array. (muM)
%   dD     - Change in [Hb] concentration in a time x 1 array. (muM)
%   dT     - Change in [HbT] concentration in a time x 1 array. (muM)
%   dmua   - Change in Absorption in a time x wavelength array. (1/cm)
%   dPh    - Change in phase time x wavelength array. (rad)
%   Ph0    - Baseline phase 1 x wavelength array. (rad)

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
    spectra=load('ext_dpf.mat');
    Oext=interp1(spectra.lambda, spectra.Oext, lambda); %1/(mM cm)
    Dext=interp1(spectra.lambda, spectra.Dext, lambda); %1/(mM cm)
    
    X=linsolve([Oext', Dext'], dmua');
    dO=X(1,:)'*1000; %uM
    dD=X(2,:)'*1000; %uM
    dT=dO+dD; %uM

end

