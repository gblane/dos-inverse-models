function [dO, dD, dT, dmua] = I2Blood_SSIn(dis, II, mua0Lam, musp0Lam, lambda, fmod, nin)
% I2Blood_SSIn Calculate hemoglobin changes from single-slope intensity measurements.
%
% [dO, dD, dT, dmua] = I2Blood_SSIn(dis, II, mua0Lam, musp0Lam, lambda, fmod, nin)
%
% Written by Giles Blaney, Ph.D. Spring 2019
%
% Inputs:
%   dis      - Source-detector distances [cm]
%   II       - Intensity data (time x distance for each wavelength) [mixed]
%   mua0Lam  - Baseline absorption coefficients [1/cm]
%   musp0Lam - Baseline reduced scattering coefficients [1/cm]
%   lambda   - Optical wavelengths [nm]
%   fmod     - (Optional, default=1.40625e8 Hz) Modulation frequency [Hz]
%   nin      - (Optional, default=1.4) Internal index of refraction [-]
%
% Outputs:
%   dO       - Change in [HbO2] concentration [muM]
%   dD       - Change in [Hb] concentration [muM]
%   dT       - Change in [HbT] concentration [muM]
%   dmua     - Change in absorption coefficient [1/cm]

    %% Setup
    if nargin<=5
        fmod=1.40625e8; %Hz
        nin=1.4;
    end
    c=2.99792458e10; %cm/sec
    omega=2*pi*fmod; %rad/sec
    nu=c/nin;
    disAvg=mean(dis);
    
    %% Params -> <L>
    Lac=ones(length(lambda), length(dis))*NaN;
    for lamInd=1:length(mua0Lam)
        mua0=mua0Lam(lamInd); %1/cm
        musp0=musp0Lam(lamInd); %1/cm
        
        for disInd=1:length(dis)
            L=complexPathLen(dis(disInd), mua0, musp0, fmod, nin);
            Lac(lamInd, disInd)=real(L); %cm
        end
    end

    %% I -> dS
    IIcal=calI_first(dis, II, mua0Lam, musp0Lam, omega, nu);

    tau=(3/2)*(mua0Lam.*musp0Lam+...
        sqrt((mua0Lam.*musp0Lam).^2+(musp0Lam.*omega/nu).^2));
    S0=-sqrt(tau);

    S=zeros(size(IIcal{1}, 1), length(lambda));
    for i=1:length(IIcal)
        lnr2II = log(IIcal{i}.*(dis.^2));
        if length(dis)==2
            S(:, i)=(lnr2II(:, 2)-lnr2II(:, 1))/(dis(2)-dis(1));
        else           
            for j=1:size(lnr2I, 1)
                pTemp=polyfit(dis,lnr2II(j, :), 1);
                S(j, i)=pTemp(1);
            end
        end
    end
    dS=S-S0;
    
    fact=-length(dis)*var(dis, 1)./sum((dis-disAvg).*Lac, 2);
    dmua=dS.*fact';

    %% dmua -> dBlood
    spectra=load('ext_dpf.mat');
    Oext=interp1(spectra.lambda, spectra.Oext, lambda); %1/(mM cm)
    Dext=interp1(spectra.lambda, spectra.Dext, lambda); %1/(mM cm)

    X=linsolve([Oext', Dext'], dmua');
    dO=X(1,:)'*1000; %uM
    dD=X(2,:)'*1000; %uM
    dT=dO+dD; %uM
    
end

function [IIcal] = calI_first(dis, II, mua, musp, omega, nu)
% [IIcal] = calI_first(dis, II, mua, musp)

    IIcal=cell(size(II));
    for i=1:length(II)
        calFact=calI(dis, mean(II{i}(1:500, :), 1), mua(i), musp(i), omega, nu);
        IIcal{i}=II{i}.*calFact;
    end

end

function [calFact] = calI(dis, I, mua, musp, omega, nu)
% [calFact] = calI(dis, I, mua, musp)

    if size(dis, 1)~=1
        dis=dis';
    end

    if size(I, 2)~=length(dis)
        I=I';
    end

    tau=(3/2)*(mua*musp+sqrt((mua*musp)^2+(musp*omega/nu)^2));
    Sl=-sqrt(tau);
    lnr2Itrue = dis*Sl;
    Itrue = exp(lnr2Itrue)./(dis.^2);
    calFact = Itrue./I;
end