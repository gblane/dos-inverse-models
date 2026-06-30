function [mua, musp, a, b] = assumeOptProp(lambda)
% assumeOptProp Estimate default tissue optical properties by wavelength.
%
% [mua, musp, a, b] = assumeOptProp(lambda)
%
% Written by Giles Blaney (DOIT-Toolbox version; Ph.D. awarded May 2022)
% Moved into dos-inverse-models and adapted to use makeE data (2026).
%
% Inputs:
%   lambda - Wavelengths [nm]. Default: [690, 830].
%
% Outputs:
%   mua    - Assumed absorption coefficient [1/mm].
%   musp   - Assumed reduced scattering coefficient [1/mm].
%   a      - Reduced scattering amplitude at 690 nm [1/mm].
%   b      - Reduced scattering power-law exponent [-].
%
% Notes:
%   The absorption estimate assumes 70% water content and follows the
%   nominal optical properties used in the DOIT-Toolbox helper:
%   mua=[0.013, 0.011] 1/mm and musp=[0.77, 0.66] 1/mm at [690, 830] nm.
%   The original helper cited doi:10.1117/1.JBO.17.8.081406.
%
% Dependencies:
%   makeE is provided by this repository and uses data in ../../data.

    if nargin < 1 || isempty(lambda)
        lambda = [690, 830];
    end
    if size(lambda, 1) == 1
        lambda = lambda.';
    end

    mua0 = [0.013, 0.011]; % 1/mm
    musp0 = [0.77, 0.66]; % 1/mm
    lambda0 = [690, 830]; % nm
    Cw = 0.7;

    E0 = makeE('ODW', lambda0);
    Oext0 = E0(:, 1).';
    Dext0 = E0(:, 2).';
    Wext0 = E0(:, 3).';

    X = (Wext0(1)*Dext0(2) - Wext0(2)*Dext0(1))/...
        (Oext0(1)*Dext0(2) - Oext0(2)*Dext0(1));
    Y = (Wext0(1)*Oext0(2) - Wext0(2)*Oext0(1))/...
        (Oext0(1)*Dext0(2) - Oext0(2)*Dext0(1));

    bloodConc = linsolve([Oext0.', Dext0.'], mua0.');
    O0 = bloodConc(1) - X*Cw;
    D0 = bloodConc(2) + Y*Cw;

    E = makeE('ODW', lambda);
    mua = O0*E(:, 1) + D0*E(:, 2) + Cw*E(:, 3);

    a = musp0(1);
    lambdaRef = lambda0(1);
    b = -log(musp0(2)/a)/log(lambda0(2)/lambdaRef);
    musp = a*(lambda/lambdaRef).^-b;
end
