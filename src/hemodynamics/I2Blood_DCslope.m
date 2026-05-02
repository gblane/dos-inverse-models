function [dO, dD, dT, dmua] = ...
    I2Blood_DCslope(dis, II, mua0, musp0, lambda)
% Giles Blaney Spring 2019
% [dO, dD, dT, dmua] = I2Blood_Islope(II, mua0, musp0, lambda, dis)
% Inputs:
%   dis    - Distances in a 1 x distnace array. (cm)
%   II     - Intensity data in a 1 x wavelength cell. Cell elements contain
%            time x distnace data.
%   mua0   - Baseline absorption in a 1 x wavelength array. (1/cm)
%   musp0  - Baseline scattering in a 1 x wavelength array. (1/cm)
%   lambda - Wavelengths in a 1 x wavelength array. (nm)
% Outputs:
%   dO     - Change in [HbO2] concentration in a time x 1 array. (muM)
%   dD     - Change in [Hb] concentration in a time x 1 array. (muM)
%   dT     - Change in [HbT] concentration in a time x 1 array. (muM)
%   dmua   - Change in Absorption in a time x wavelength array. (1/cm)

    IIcal=calI_first(dis, II, mua0, musp0);

    S=zeros(size(IIcal{1}, 1), length(lambda));
    for i=1:length(IIcal)
        lnr2I=log(dis.^2.*IIcal{i});
        if length(dis)==2
            S(:,i) = (lnr2I(:,1)-lnr2I(:,2))./(dis(1)-dis(2));
        else           
            for j=1:size(lnr2I, 1)
                pTemp=polyfit(dis, lnr2I(j, :), 1);
                S(j, i)=pTemp(1);
            end
        end        
    end
    
    mua=S.^2./(3*musp0);
    dmua=mua-mua0;
    
    spectra=load('ext_dpf.mat');
    Oext=interp1(spectra.lambda, spectra.Oext, lambda); %1/(mM cm)
    Dext=interp1(spectra.lambda, spectra.Dext, lambda); %1/(mM cm)
    
    X=linsolve([Oext', Dext'], dmua');
    dO=X(1,:)'*1000; %uM
    dD=X(2,:)'*1000; %uM
    dT=dO+dD; %uM
    
end

function [IIcal] = calI_first(dis, II, mua, musp)
% [IIcal] = calI_first(dis, II, mua, musp)

    IIcal=cell(size(II));
    for i=1:length(II)
        calFact=calI(dis, mean(II{i}(1:10, :), 1), mua(i), musp(i));
        IIcal{i}=II{i}.*calFact;
    end

end

function [calFact] = calI(dis, I, mua, musp)
% [calFact] = calI(dis, I, mua, musp)

    if size(dis, 1)~=1
        dis=dis';
    end

    if size(I, 2)~=length(dis)
        I=I';
    end

    Sl = -sqrt(3*mua*musp);
    lnr2Itrue = polyval([Sl 0],dis);
    Itrue = exp(lnr2Itrue)./(dis.^2);
    calFact = Itrue./I;
end