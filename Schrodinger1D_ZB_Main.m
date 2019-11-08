%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 8November2019, lne %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves the Schrodinger equation with m(E,z) using FDM (Finite Difference Method) 
%
% The non-parabolicity is implemented via the kp model using the k.p 2 bands Kane model
% for the conduction band while it is using the k.p 6bands Luttinger-Kohn model for
% the valence bands
% The HH band remains ALL THE TIME parabolic in ZB-001 growth and even with strain
% A strain model is included. It basically shifts the conduction and valence band edge
% but also influences the coupling bw LH/SO in the kp 6bands.
% The strain is mainly interesting for InGaAs/GaAs heterostructures

% -> II-VI and cubic nitride material parameters are available but should
% be grabt in the "Library.m" file
% -> Additionnal material can be added in the "materialDB_ZB.csv" file
% -> Additionnal Ternary alloy material can be added in the "TernaryAlloy.m" file
% -> Additionnal Quaternary alloy material can be added in the "QuaternaryAlloy.m" file

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;         % Planck constant J.s
hbar=h/(2*pi);
e=1.602176487E-19;        % charge de l electron Coulomb
m0=9.10938188E-31;        % electron mass kg
c=2.99792458e8;           % speed of light (m/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StrainModel          = 1; % Activate Strain model
Display_IB_Results   = 1; % Switch to print or not the IB  dipoles on the shell
Display_ISB_Results  = 1; % Switch to print or not the ISB dipoles on the shell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=3;                      % number of solution asked per model
ScF=0.1;                  % scaling factor to plot the wave function [Without Dimension]
dz=1E-10;                 % resolution of the grid [m]
F0=0;%1e7;                % Electric field [Volt/meter]
T=300;                    % Temperature [Kelvin], react on the band gap Eg only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library;                  % load material parameter DB from "materialDB_ZB.csv"
ExtractParameters;        % extract parameter from the Library
TernaryAlloy;             % compute the ternary alloy
QuaternaryAlloy;          % compute the quaternary alloy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% import the layer structure file %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE !!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt   = M(:,end)*1e-9;    % conversion of the length from Angstrom to meter

Egt  = M(:,idx_Eg6c) - (M(:,idx_alphaG)*T^2) ./ (T+M(:,idx_betaG));   %Eg = Eg0 - (a*T.^2)./(T + b);
VBOt = M(:,idx_VBO);
CBOt = Egt+VBOt;         % CBO form band gap difference and temperature
Dsot = M(:,idx_Dso);     % Spin-Orbit shift band parameter
%Ft   = M(:,idx_F);      % Gammac Luttinger parameter for the electron (used for k.p 8bands only)
g1t  = M(:,idx_g1);      % Gamma1 Luttinger parameter
g2t  = M(:,idx_g2);      % Gamma2 Luttinger parameter
g3t  = M(:,idx_g3);      % Gamma3 Luttinger parameter

EPt_K= M(:,idx_EP_K);    % EP Kane
%EPt_L= M(:,idx_EP_L);   % EP Luttinger (used for k.p 8bands only)
%Epsit= M(:,idx_Epsi);   %(used for Poisson solver only)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Strain Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

at  = M(:,idx_a);           % lattice parameter
act = M(:,idx_ac);          % Conduction band strain offset parameter
avt = M(:,idx_av);          % Valence band strain offset parameter
bvt = M(:,idx_bv);          % Valence band strain offset parameter
c11t = M(:,idx_c11);        % strain parameter
c12t = M(:,idx_c12);        % strain parameter

a0   = substrate(idx_a);

if StrainModel == 1
  exxt =  (a0-at)/a0; % eyyt =  exxt;
  ezzt = -2*c12t./c11t.*exxt;
else
  exxt =  (a0-at)/a0 * 0; % eyyt =  exxt;
  ezzt = -2*c12t./c11t.*exxt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z and a lot of other parameters

z=0; V0=CBOt(1); Eg=Egt(1);Dso=Dsot(1);
%Epsi=Epsit(1); F=Ft(1);
g1=g1t(1); g2=g2t(1); g3=g3t(1); EP_K=EPt_K(1); %EP_L=EPt_L(1); %me(1)=met(1);% mhh(1)=mhh_t(1);
ac=0; av=0; bv=0; exx=0; ezz=0;

for i=1:length(zt)
    zv  =  (z(end)+dz) : dz : (z(end) + zt(i)) ;
    z   = [ z  zv ];
    V0  = [ V0     ones(size(zv)) * CBOt(i)  ];
    Eg  = [ Eg     ones(size(zv)) * Egt(i)   ];
    EP_K= [ EP_K   ones(size(zv)) * EPt_K(i) ];
%   EP_L= [ EP_L   ones(size(zv)) * EPt_L(i) ];
    Dso = [ Dso    ones(size(zv)) * Dsot(i)  ];
%   F   = [ F      ones(size(zv)) * Ft(i)    ];
    g1  = [ g1     ones(size(zv)) * g1t(i)   ];
    g2  = [ g2     ones(size(zv)) * g2t(i)   ];
    g3  = [ g3     ones(size(zv)) * g3t(i)   ];
    ac  = [ ac     ones(size(zv)) * act(i)   ];
    av  = [ av     ones(size(zv)) * avt(i)   ];
    bv  = [ bv     ones(size(zv)) * bvt(i)   ];
    exx = [ exx    ones(size(zv)) * exxt(i)  ];
    ezz = [ ezz    ones(size(zv)) * ezzt(i)  ];
end

V0=V0-min(V0);             % Shift the band in order to get the bottom of the well at zero
V0=(F0*z)+V0;              % adding the electric field to the potential

eyy = exx;
DCBO   = -abs(ac).*(exx+eyy+ezz) ;                      % shift of the CB due to strain
DVBOHH = +abs(av).*(exx+eyy+ezz) - abs(bv).*(exx-ezz) ; % shift of the VB-HH due to strain
DVBOLH = +abs(av).*(exx+eyy+ezz) + abs(bv).*(exx-ezz) ; % shift of the VB-LH due to strain
DVBOSO = +abs(av).*(exx+eyy+ezz) ;                      % shift of the VB-SO due to strain
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ec,psic] = Schrod_2bands_Kane_f(z,V0,Eg,EP_K,Dso,n,ac,av,bv,exx,ezz);

% m=m(z), take care, the mhh in the z-direction is different from mhh in 3D
mhhL= 1 ./ (g1-2*g2);
[Ehh,psihh] = Schrod_1band_f(z,-(V0-Eg+DVBOHH),mhhL,n);  % m = m(z) the HH are all the time parabolic in ZB-001, even with strain!
Ehh=-Ehh;
[Elh,psilh] = Schrod_2bands_Luttinger_Kohn_f(z,V0,Eg,Dso,g1,g2,g3,n,av,bv,exx,ezz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale_PSI;
computesISBdipoles;
computesIBoverlap;

if Display_IB_Results == 1
  PrintIBResults;
end
if Display_ISB_Results == 1
  PrintISBResults;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[10 100 1000 700]);
subplot(1,1,1,'fontsize',15)
hold on;grid on;

xscale=[z(1) z(end)]*1e9;
shift=0;
yscale=[min(V0-Eg)-0.4 max(V0)+0.1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(z*1e9,V0        + shift , 'b-','linewidth',2)
plot(z*1e9,V0-Eg     + shift , 'b-' ,'linewidth',2)
plot(z*1e9,V0-Eg-Dso + shift , 'b-','linewidth',2)
  
if StrainModel == 1
  plot(z*1e9,V0   +DCBO       + shift , 'b--','linewidth',2)
  plot(z*1e9,V0-Eg+DVBOHH     + shift , 'm--' ,'linewidth',2)
  plot(z*1e9,V0-Eg+DVBOLH     + shift , 'c--' ,'linewidth',2)
  plot(z*1e9,V0-Eg-Dso+DVBOSO + shift , 'c--' ,'linewidth',2)
end


for i=1:length(Ec)
    plot(z*1e9,PSIc(:,i)+ shift,'color','r','linewidth',2)
end

for i=1:length(Ehh)
  plot(z*1e9,PSIhh(:,i)+ shift,'color','m','linewidth',2)
end
for i=1:length(Elh)
  plot(z*1e9,PSIlh(:,i)+ shift,'color','c','linewidth',2)
end

xlabel('z (nm)');
ylabel('Energy (eV)');

if StrainModel == 1
  title(strcat('T=',num2str(T),'K ; with STRAIN'))
else
  title(strcat('T=',num2str(T),'K ; without STRAIN'))
end

xlim(xscale)
ylim(yscale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%