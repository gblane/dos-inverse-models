function [f, J, mod, fram] = fun_marq2DE_six_param(x0, fpar, data, flags, en)
% fun_marq2DE_six_param Calculate residuals and Jacobian for two-layer Marquardt recovery.
%
% [f, J, mod, fram] = fun_marq2DE_six_param(x0, fpar, data, flags, en)
%
% Written by Giles Blaney, Ph.D. (Originally by Bertan H. May 2015)
%
% Inputs:
%   x0    - Fitting parameters [mua1, musp1, L, mua2, musp2, af] [mixed]
%   fpar  - Fitting parameters indices [-]
%   data  - Experimental data [mm, rad, arb, rad, arb] [mixed]
%   flags - Vector of booleans controlling parameters to fit [-]
%   en    - Zeroth order Bessel function roots [-]
%
% Outputs:
%   f     - Vector of residuals [unitless]
%   J     - Jacobian matrix [mixed]
%   mod   - Forward model output [mixed]
%   fram  - Frame information [-]
format long
r0 = data(:,1); ph0 = data(:,2); ac0 = data(:,3); phs = data(:,4); acs = data(:,5);
% Units are in mm-1 and mm
ua1 = x0(1); us1 = x0(2); L0 = x0(3);
ua2 = x0(4); us2 = x0(5); af = x0(6);
mu = [ua1,ua2,us1,us2];
% [AC, phase] = Reflectance(mu,L0,af,r0);
[AC, phase] = Reflectance_Ceph(mu,L0,r0,en);
%
%APRIL 15: I HAVE CHANGED 'REFLECTANCE' WITH 'REFLECTANCE_Ceph'
%%%% Numerical Derivatives %%%%
derVec = [1e-04,1e-02,0.001,1e-04,1e-02,1e-04].*flags;
dua1 = derVec(1); dus1 = derVec(2); 
dua2 = derVec(4); dus2 = derVec(5);
daf  = derVec(6); dL0  = derVec(3);
% dua1 = 1e-04; dus = 1e-02; daf = 1e-04; dL0 = 0.001;
% 1) Perturbed Layer 1 ua
mu = [ua1+dua1,ua2,us1,us2];
[AC_ua1p, phase_ua1p] = Reflectance_Ceph(mu,L0,r0,en);
mu = [ua1-dua1,ua2,us1,us2];
[AC_ua1m, phase_ua1m] = Reflectance_Ceph(mu,L0,r0,en);
% 2) Perturbed Layer 1 us'
mu = [ua1,ua2,us1+dus1,us2];
[AC_us1p, phase_us1p] = Reflectance_Ceph(mu,L0,r0,en);
mu = [ua1,ua2,us1-dus1,us2];
[AC_us1m, phase_us1m] = Reflectance_Ceph(mu,L0,r0,en);
% 3) Perturbed Layer 1 Thickness
mu = [ua1,ua2,us1,us2];
[AC_L0p, phase_L0p] = Reflectance_Ceph(mu,L0+dL0,r0,en);
[AC_L0m, phase_L0m] = Reflectance_Ceph(mu,L0-dL0,r0,en);
% 4) Perturbed Layer 2 ua
mu = [ua1,ua2+dua2,us1,us2];
[AC_ua2p, phase_ua2p] = Reflectance_Ceph(mu,L0,r0,en);
mu = [ua1,ua2-dua2,us1,us2];
[AC_ua2m, phase_ua2m] = Reflectance_Ceph(mu,L0,r0,en);
% 5) Perturbed Layer 2 us'
mu = [ua1,ua2,us1,us2+dus2];
[AC_us2p, phase_us2p] = Reflectance_Ceph(mu,L0,r0,en);
mu = [ua1,ua2,us1,us2-dus2];
[AC_us2m, phase_us2m] = Reflectance_Ceph(mu,L0,r0,en);
% 6) Perturbed Amplitude Factor
mu = [ua1,ua2,us1,us2];
[AC_afp, phase_afp] = Reflectance_Ceph(mu,L0,r0,en);
[AC_afm, phase_afm] = Reflectance_Ceph(mu,L0,r0,en);

der_acmod(:,1) = (1./acs).*((AC_ua1p-AC_ua1m)/(2*dua1));       
der_acmod(:,2) = (1./acs).*((AC_us1p-AC_us1m)/(2*dus1));
der_acmod(:,3) = (1./acs).*((AC_L0p-AC_L0m)/(2*dL0));  
der_acmod(:,4) = (1./acs).*((AC_ua2p-AC_ua2m)/(2*dua2));       
der_acmod(:,5) = (1./acs).*((AC_us2p-AC_us2m)/(2*dus2));     
der_acmod(:,6) = (1./acs).*((AC_afp-AC_afm)/(2*daf));         
der_phmod(:,1) = (1./phs).*((phase_ua1p-phase_ua1m)/(2*dua1));
der_phmod(:,2) = (1./phs).*((phase_us1p-phase_us1m)/(2*dus2));
der_phmod(:,3) = (1./phs).*((phase_L0p-phase_L0m)/(2*dL0));    
der_phmod(:,4) = (1./phs).*((phase_ua2p-phase_ua2m)/(2*dua2));
der_phmod(:,5) = (1./phs).*((phase_us2p-phase_us2m)/(2*dus2));
der_phmod(:,6) = (1./phs).*((phase_afp-phase_afm)/(2*daf));

% [AC phase]
costac = (AC-ac0)./acs; costac = costac'*costac;
costph = (phase-ph0)./phs; costph = costph'*costph;
er = [costac costph];

f = [(AC-ac0)./acs ; (phase-ph0)./phs]; %Cost vector
J = [der_acmod;der_phmod];              %Jacobian matrix of the cost
J(isnan(J)) = 0;
% W = diag([1./(acs);1./(phs)]);

% A = J'*J;
% diag(A)

%%%% Video to keep track of progress %%%%
mod = [r0,phase,AC];
%%%comment the next lines to the end for not plotting fits

% v = figure(1); clf;
% set(gcf,'Units','pixels','Position', [300 270 800 300]);
% % set(gca,'Color','None');
% set(gca,'FontSize',12)
% subplot(121), semilogy(r0,ac0,'bo'), hold on, semilogy(r0,AC,'r'); 
% subplot(122), plot(r0,ph0,'bo'), hold on, plot(r0,phase,'r'), ylim([0 1.2])
% fram = getframe(v);

fram=[];
