function [LY] = calcPathLen_datTyp(R, L, Rcw, Lcw, datTyp)
% calcPathLen_datTyp Calculate the pathlength of a specific data type.
%
% [LY] = calcPathLen_datTyp(R, L, Rcw, Lcw, datTyp)
%
% Written by Giles Blaney (Summer 2020; Ph.D. awarded May 2022) (Updated Winter 2023)
%
% Inputs:
%   R      - Complex reflectance [1/mm^2]
%   L      - Complex pathlength [mm]
%   Rcw    - Real CW reflectance [1/mm^2]
%   Lcw    - Real CW pathlength [mm]
%   datTyp - String for data type for which to calculate pathlength [-]
%
% Outputs:
%   LY     - Pathlength of 'datTyp' [mm]

    switch datTyp
        case 'I' % ln(r^2 |R|)
            LY=real(L);
            
        case 'P' % angle(R)
            LY=imag(L);

        case 'C' % ln(r^2 I)
            LY=Lcw;
            
        case 'Re' % ln(r^2 Re(R))
            LY=real(L)-(imag(R)./real(R)).*imag(L);
            
        case 'Im' % ln(r^2 Im(R))
            LY=real(L)+(real(R)./imag(R)).*imag(L);
            
        case 'ReN' % Re(R/Rcw);
            L_Re=calcPathLen_datTyp(R, L, Rcw, Lcw, 'Re');
            LY=real(R./Rcw).*(L_Re-Lcw);

        case 'OmReN' % 1-Re(R/Rcw)
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            LY=-L_ReN;
            
        case 'ImN' % Im(R/Rcw)
            L_Im=calcPathLen_datTyp(R, L, Rcw, Lcw, 'Im');
            LY=imag(R./Rcw).*(L_Im-Lcw);
        
        case {'ReNpImN', 'ReNpImNmO'} % Re(R/Rcw)+Im(R/Rcw) & Re(R/DC)+Im(R/DC)-1
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            LY=L_ReN+L_ImN;







        case {'ImNmReN', 'ImNmReNpO'} % Im(R/Rcw)-Re(R/Rcw) & Im(R/DC)-Re(R/DC)+1
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            LY=L_ImN-L_ReN;
            
            


            


        case {'ReNpP', 'ReNpPmO'} % Re(R/Rcw)+phi & Re(R/Rcw)+phi-1
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ReN+L_P;
            






        case 'ReNmP' % Re(R/Rcw)-phi
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ReN-L_P;
            warning('Deprecated data-type');
            
        case {'PmReN', 'PmReNpO'} % phi-Re(R/Rcw) & phi-Re(R/Rcw)+1
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_P-L_ReN;
            






        case 'ImNpP' % Im(R/Rcw)+phi
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ImN+L_P;
            
        case 'ImNmP' % Im(R/Rcw)-phi
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ImN-L_P;
            warning('Deprecated data-type');

        case 'PmImN' % phi-Im(R/Rcw)
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            LY=L_P-L_ImN;

        case 'ReNpP2' % Re(R/Rcw)+phi^2/2
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ReN+angle(R).*L_P;
            warning('Deprecated data-type');
            
        case 'OmReNpP2' % 1-[Re(R/Rcw)+phi^2/2]
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=-(L_ReN+angle(R).*L_P);
            
        case 'IoC' % I/C
            LY=abs(R)./Rcw.*(real(L)-Lcw);

        case 'OmIoC' % 1-I/C
            LY=-abs(R)./Rcw.*(real(L)-Lcw);

        case 'OmImNoP' %1-Im(R/Rcw)/phi
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            LY=-1./angle(R).*L_ImN+1./angle(R).^2.*imag(R./Rcw).*L_P;

        case 'mReNpImNmO' % Re(R/Rcw)+Im(R/Rcw) & Re(R/DC)+Im(R/DC)-1
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            LY=-L_ReN-L_ImN;
            
        otherwise
            error('Unknown data type');
    end
end