function [E]=makeE(chroms, lambda)
% [E]=makeE(chroms, lambda)
% Giles Blaney Spring 2021
% 
% Inputs:   - chroms: String of chromophores to include in E. 
%                     Available chromophores:
%                     - O: Oxyhemoglobin
%                     - D: Deoxyhemoglobin
%                     - W: Water
%                     - L: Lipid
%                     - C: Collagen
%                     - CCOo: Oxidized cytochrome c oxidase
%                     - CCOr: reduced cytochrome c oxidase
%                     (Default: 'OD')
%                     Spaces are ignored.
%           - lambda: Vectors of wavelengths (nm).
%                     (Default: [830, 690])
% 
% Output:   - E: Extinction coefficient matrix.
%                Units: 1/(mm uM) for O, D, CCOo, and CCOr
%                       1/mm for W, L, and C
%                size(E)=[length(lambda), length(chroms)];
%                Defined as mua=E*C
%                Order of C is defined by order in chroms input
    
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
                blood=load('Bext.mat');
                Oext=interp1(blood.Blambda, blood.Oext, lambda,...
                    'linear', 'extrap');
                Dext=interp1(blood.Blambda, blood.Dext, lambda,...
                    'linear', 'extrap');
                Oext(Oext<0)=0;
                Dext(Dext<0)=0;
                E=[E, Oext+Dext];
            case 'O' %Oxy
                blood=load('Bext.mat');
                Oext=interp1(blood.Blambda, blood.Oext, lambda,...
                    'linear', 'extrap');
                Oext(Oext<0)=0;
                E=[E, Oext];
            case 'D' %Deoxy
                blood=load('Bext.mat');
                Dext=interp1(blood.Blambda, blood.Dext, lambda,...
                    'linear', 'extrap');
                Dext(Dext<0)=0;
                E=[E, Dext];
            case 'W' %Water
                water=load('Wext.mat');
                Wext=interp1(water.Wlambda, water.Wmua, lambda,...
                    'linear', 'extrap');
                Wext(Wext<0)=0;
                E=[E, Wext];
            case 'L' %Lipid
                lipid=load('Lext.mat');
                Lext=interp1(lipid.Llambda, lipid.Lmua, lambda,...
                    'linear', 'extrap');
                Lext(Lext<0)=0;
                E=[E, Lext];
            case 'C' %Collagen or CCO
                if strcmp(chroms(1:3), 'CCO') %CCO
                    cco=load('CCOext.mat');
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
                    chroms(1:3)=[];
                else %Collagen
                    col=load('Cext.mat');
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