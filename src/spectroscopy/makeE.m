function [E] = makeE(chroms, lambda)
% makeE Generate the extinction coefficient matrix for specified chromophores.
%
% [E] = makeE(chroms, lambda)
%
% Written by Giles Blaney (Spring 2021; Ph.D. awarded May 2022)
%
% Inputs:
%   chroms - String of chromophores to include (e.g., 'ODWL') [-]
%   lambda - Vector of wavelengths [nm]
%
% Outputs:
%   E      - Extinction coefficient matrix [1/(mm uM) or 1/mm]
    
    if nargin<=0
        chroms='OD';
    end
    if nargin<=1
        lambda=[830, 690];
    end
    if size(lambda, 1)==1
        lambda=lambda';
    end
    
    chroms=chroms(~isspace(chroms));
    
    E=[];
    while ~isempty(chroms)
        switch chroms(1)
            case 'T' %Total
                blood=loadSpectra('Bext.mat');
                Oext=interp1(blood.Blambda, blood.Oext, lambda,...
                    'linear', 'extrap');
                Dext=interp1(blood.Blambda, blood.Dext, lambda,...
                    'linear', 'extrap');
                Oext(Oext<0)=0;
                Dext(Dext<0)=0;
                E=[E, Oext+Dext];
            case 'O' %Oxy
                blood=loadSpectra('Bext.mat');
                Oext=interp1(blood.Blambda, blood.Oext, lambda,...
                    'linear', 'extrap');
                Oext(Oext<0)=0;
                E=[E, Oext];
            case 'D' %Deoxy
                blood=loadSpectra('Bext.mat');
                Dext=interp1(blood.Blambda, blood.Dext, lambda,...
                    'linear', 'extrap');
                Dext(Dext<0)=0;
                E=[E, Dext];
            case 'W' %Water
                water=loadSpectra('Wext.mat');
                Wext=interp1(water.Wlambda, water.Wmua, lambda,...
                    'linear', 'extrap');
                Wext(Wext<0)=0;
                E=[E, Wext];
            case 'L' %Lipid
                lipid=loadSpectra('Lext.mat');
                Lext=interp1(lipid.Llambda, lipid.Lmua, lambda,...
                    'linear', 'extrap');
                Lext(Lext<0)=0;
                E=[E, Lext];
            case 'C' %Collagen or CCO
                if length(chroms)>=4 && strncmp(chroms, 'CCO', 3) %CCO
                    cco=loadSpectra('CCOext.mat');
                    switch chroms(4)
                        case 'o' %CCOo
                            oCCOext=...
                                interp1(cco.CCOlambda, cco.oCCOext, lambda,...
                                'linear', 'extrap');
                            oCCOext(oCCOext<0)=0;
                            E=[E, oCCOext];
                        case 'r' %CCOr
                            rCCOext=...
                                interp1(cco.CCOlambda, cco.rCCOext, lambda,...
                                'linear', 'extrap');
                            rCCOext(rCCOext<0)=0;
                            E=[E, rCCOext];
                        otherwise
                            warning('Unknown state of CCO, ignored');
                    end
                    chroms=chroms(4:end);
                else %Collagen
                    col=loadSpectra('Cext.mat');
                    Cext=interp1(col.Clambda, col.Cmua, lambda,...
                        'linear', 'extrap');
                    Cext(Cext<0)=0;
                    E=[E, Cext];
                end
            otherwise
                warning(['Unknown chromophore' chroms(1) ', ignored']);
        end
        chroms(1)=[];
    end
end

function spectra = loadSpectra(fileName)
% loadSpectra Load spectroscopy data from this repository's data folder.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    spectra = load(fullfile(repoRoot, 'data', fileName));
end
