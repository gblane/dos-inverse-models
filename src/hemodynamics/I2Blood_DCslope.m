function [dO, dD, dT, dmua] = I2Blood_DCslope(dis, II, mua0, musp0, lambda)
% I2Blood_DCslope Calculate hemoglobin changes from DC slope measurements.
%
% [dO, dD, dT, dmua] = I2Blood_DCslope(dis, II, mua0, musp0, lambda)
%
% Written by Giles Blaney, Ph.D. Spring 2019
%
% Inputs:
%   dis    - Source-detector distances [cm]
%   II     - Intensity data (time x distance for each wavelength) [mixed]
%   mua0   - Baseline absorption coefficients [1/cm]
%   musp0  - Baseline reduced scattering coefficients [1/cm]
%   lambda - Optical wavelengths [nm]
%
% Outputs:
%   dO     - Change in [HbO2] concentration [muM]
%   dD     - Change in [Hb] concentration [muM]
%   dT     - Change in [HbT] concentration [muM]
%   dmua   - Change in absorption coefficient [1/cm]

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