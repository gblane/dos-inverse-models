function [dO, dD, dT, dmua] = ...
    I2Blood_SSPh(dis, PhPh, mua0Lam, musp0Lam, lambda, fmod, nin)
% Giles Blaney Spring 2019
% [dO, dD, dT, dmua] = I2Blood_SSPh(II, mua0, musp0, lambda, dis)
% Inputs:
%   dis    - Distances in a 1 x distnace array. (cm)
%   PhPh   - Phase data in a 1 x wavelength cell. Cell elements contain
%            time x distnace data.
%   mua0   - Baseline absorption in a 1 x wavelength array. (1/cm)
%   musp0  - Baseline scattering in a 1 x wavelength array. (1/cm)
%   lambda - Wavelengths in a 1 x wavelength array. (nm)
%   fmod   - (OPTIONAL, default=1.40625e8 Hz) Modulation frequecy. (Hz)
%   nin    - (OPTIONAL, default=1.4) Internal index of refraction. (-)
% Outputs:
%   dO     - Change in [HbO2] concentration in a time x 1 array. (muM)
%   dD     - Change in [Hb] concentration in a time x 1 array. (muM)
%   dT     - Change in [HbT] concentration in a time x 1 array. (muM)
%   dmua   - Change in Absorption in a time x wavelength array. (1/cm)

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
    Lph=ones(length(lambda), length(dis))*NaN;
    for lamInd=1:length(mua0Lam)
        mua0=mua0Lam(lamInd); %1/cm
        musp0=musp0Lam(lamInd); %1/cm
        
        for disInd=1:length(dis)
            L=complexPathLen(dis(disInd), mua0, musp0, fmod, nin);
            Lph(lamInd, disInd)=imag(L); %cm
        end
    end

    %% Ph -> dS
    PhPhcal=calPh_first(dis, PhPh, mua0Lam, musp0Lam, omega, nu);

    tau=(3/2)*(mua0Lam.*musp0Lam+sqrt((mua0Lam.*musp0Lam).^2+(musp0Lam.*omega/nu).^2));
    S0=(3/2).*musp0Lam.*(omega/nu).*(1./sqrt(tau));

    S=zeros(size(PhPhcal{1}, 1), length(lambda));
    for i=1:length(PhPhcal)
        if length(dis)==2
            S(:, i)=(PhPhcal{i}(:, 2)-PhPhcal{i}(:, 1))/(dis(2)-dis(1));
        else           
            for j=1:size(lnr2I, 1)
                pTemp=polyfit(dis, PhPhcal{i}(j, :), 1);
                S(j, i)=pTemp(1);
            end
        end
    end
    dS=S-S0;
    
    fact=-length(dis)*var(dis, 1)./sum((dis-disAvg).*Lph, 2);
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

%% Functions
function [PhPhcal] = calPh_first(dis, PhPh, mua, musp, omega, nu)

    PhPhcal=cell(size(PhPh));
    for i=1:length(PhPh)
        cal=calPh(dis, wrapTo2Pi(circ_mean(PhPh{i}(1:10, :))),...
            mua(i), musp(i), omega, nu);
        PhPhcal{i}=PhPh{i}-cal;
    end

end

function [cal] = calPh(dis, Ph, mua, musp, omega, nu)

    if size(dis, 1)~=1
        dis=dis';
    end

    if size(Ph, 2)~=length(dis)
        Ph=Ph';
    end
    
    tau=(3/2)*(mua*musp+sqrt((mua*musp)^2+(musp*omega/nu)^2));

    Sphi=(3/2)*musp*(omega/nu)*(1/sqrt(tau));
    
    Ph_exp=dis*Sphi;
    cal=Ph-Ph_exp;
end