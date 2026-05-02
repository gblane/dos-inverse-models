function [Y] = calcData_datTyp(rho, R, Rcw, datTyp)
% calcData_datTyp Calculate various data types from complex and CW reflectance.
%
% [Y] = calcData_datTyp(rho, R, Rcw, datTyp)
%
% Written by Giles Blaney, Ph.D. Summer 2020 (Updated Winter 2023)
%
% Inputs:
%   rho    - Source-detector distance [mm]
%   R      - Complex reflectance [1/mm^2]
%   Rcw    - Real CW reflectance [1/mm^2]
%   datTyp - String for data type for which to calculate data value [-]
%
% Outputs:
%   Y      - Calculated data value [mixed]

    switch datTyp
        case 'I' % ln(r^2 |R|)
            Y=log(rho.^2.*abs(R));
            
        case 'P' % angle(R)
            Y=angle(R);

        case 'C' % ln(r^2 I)
            Y=log(rho.^2.*Rcw);
            
        case 'Re' % ln(r^2 Re(R))
            Y=log(rho.^2.*real(R));
            
        case 'Im' % ln(r^2 Im(R))
            Y=log(rho.^2.*imag(R));
            
        case 'ReN' % Re(R/Rcw);
            Y=real(R./Rcw);
        

        case 'OmReN' % 1-Re(R/Rcw)
            Y=1-real(R./Rcw);

            
        case 'ImN' % Im(R/Rcw)
            Y=imag(R./Rcw);
            
        
        case 'ReNpImN' % Re(R/DC)+Im(R/DC)
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            Y=ReN+ImN;
            warning('Deprecated data-type');

        case 'ReNpImNmO' % Re(R/DC)+Im(R/DC)-1
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            Y=ReN+ImN-1;

        case 'ImNmReN' % Im(R/Rcw)-Re(R/Rcw)
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            Y=ImN-ReN;
            warning('Deprecated data-type');
            
        case 'ImNmReNpO' % Im(R/DC)-Re(R/DC)+1
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            Y=ImN-ReN+1;
            
        case 'ReNpP' % Re(R/Rcw)+phi
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ReN+P;
            warning('Deprecated data-type');

        case 'ReNpPmO' % Re(R/Rcw)+phi-1
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ReN+P-1;

        case 'ReNmP' % Re(R/Rcw)-phi
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ReN-P;
            warning('Deprecated data-type');        
            
        case 'PmReN' % phi-Re(R/Rcw);
            P=calcData_datTyp(rho, R, Rcw, 'P');
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            Y=P-ReN;

        case 'PmReNpO' % phi-Re(R/Rcw)+1
            P=calcData_datTyp(rho, R, Rcw, 'P');
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            Y=P-ReN+1;
%             warning('Deprecated data-type');
            
        case 'ImNpP' % Im(R/Rcw)+phi
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ImN+P;
            
        case 'ImNmP' % Im(R/Rcw)-phi 
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ImN-P;
            warning('Deprecated data-type');
            
        case 'PmImN' % phi-Im(R/Rcw)
            P=calcData_datTyp(rho, R, Rcw, 'P');
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            Y=P-ImN;

        case 'ReNpP2' % Re(R/Rcw)+phi^2/2
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ReN+P.^2/2;
            warning('Deprecated data-type');
            
        case 'OmReNpP2' % 1-[Re(R/Rcw)+phi^2/2]
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=1-(ReN+P.^2/2);
            
        case 'IoC' % AC/DC
            Y=abs(R./Rcw);
            
        case 'OmIoC' % 1-AC/DC
            Y=1-abs(R./Rcw);

        case 'OmImNoP' % 1-Im(R/Rcw)/phi
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=1-ImN./P;

        case 'mReNpImNmO' % 1-(Re(R/DC)+Im(R/DC))
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            Y=1-ReN-ImN;
                 
        otherwise
            error('Unknown data type');
    end
end