function [muspExtrap, b]=extrapMUSP(lambda, musp, lambdaExtrap)
% [muspExtrap, b]=extrapMUSP(lambda, musp, lambdaExtrap)
% Giles Blaney Spring 2021
% Inputs:   - lambda: Vector of wavelengths (2 long)
%           - musp: Vector of reduced scattering (2 long)
%           - lambdaExtrap: Vector of wavelengths to extrapolate along
% 
% Outputs:  - muspExtrap: Extrapolated reduced scattering
%                         (length(lambdaExtrap) long)
%           - b: Power for musp(lambda)=musp0*(lambda/lambda0).^-b
    
    b=log(musp(1)/musp(2))/log(lambda(2)/lambda(1));
    
    if b<=0 % If b is negative assume constant scattering
        b=0;
        muspExtrap=ones(size(lambdaExtrap))*mean(musp);
    else
        muspExtrap=musp(1)*(lambdaExtrap/lambda(1)).^-b;
    end
end