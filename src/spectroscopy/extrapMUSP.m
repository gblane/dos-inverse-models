function [muspExtrap, b] = extrapMUSP(lambda, musp, lambdaExtrap)
% extrapMUSP Extrapolate reduced scattering coefficients using a power law.
%
% [muspExtrap, b] = extrapMUSP(lambda, musp, lambdaExtrap)
%
% Written by Giles Blaney (Spring 2021; Ph.D. awarded May 2022)
%
% Inputs:
%   lambda       - Vector of wavelengths (2 long) [nm]
%   musp         - Vector of reduced scattering coefficients (2 long) [1/mm]
%   lambdaExtrap - Vector of wavelengths to extrapolate along [nm]
%
% Outputs:
%   muspExtrap   - Extrapolated reduced scattering coefficients [1/mm]
%   b            - Power law exponent for musp(lambda)=musp0*(lambda/lambda0).^-b [-]
    
    b=log(musp(1)/musp(2))/log(lambda(2)/lambda(1));
    
    if b<=0 % If b is negative assume constant scattering
        b=0;
        muspExtrap=ones(size(lambdaExtrap))*mean(musp);
    else
        muspExtrap=musp(1)*(lambdaExtrap/lambda(1)).^-b;
    end
end