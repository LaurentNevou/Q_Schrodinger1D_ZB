%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 8November2019, lne %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program does only a big loop of the SchrodingerSolver over one Layer thickness
% It plots various results and parameters over the Layer thickness

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
Display_IB_Results   = 0; % Switch to print or not the IB  dipoles on the shell
Display_ISB_Results  = 0; % Switch to print or not the ISB dipoles on the shell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=2;                      % number of solution asked per model
ScF=0.1;                  % scaling factor to plot the wave function [Without Dimension]
dz=1E-10;                 % resolution of the grid [m]
F0=0;%1e7;                % Electric field [Volt/meter]
T=300;                    % Temperature [Kelvin], react on the band gap Eg only
L=1:1:20;                 % Thickness of the Layer on which is done the loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library;                  % load material parameter DB from "materialDB_ZB.csv"
ExtractParameters;        % extract parameter from the Library
TernaryAlloy;             % compute the ternary alloy
QuaternaryAlloy;          % compute the quaternary alloy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=1:length(L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% import the layer structure file %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%substrate=InP;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%AlInAs  10
%InGaAs  L(l)
%AlInAs  10
%];

substrate=GaAs;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
M=[
AlGaAs30  10
GaAs  L(l)
AlGaAs30  10
];

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
ac=act(1); av=avt(1); bv=bvt(1); exx=exxt(1); ezz=ezzt(1);

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
EcL(:,l)=Ec;

% m=m(z), take care, the mhh in the z-direction is different from mhh in 3D
mhhL= 1 ./ (g1-2*g2);  
[Ehh,psihh] = Schrod_1band_f(z,-(V0-Eg+DVBOHH),mhhL,n);  % m = m(z) the HH are all the time parabolic in ZB-001, even with strain!
Ehh=-Ehh;
EhhL(:,l)=Ehh;

[Elh,psilh] = Schrod_2bands_Luttinger_Kohn_f(z,V0,Eg,Dso,g1,g2,g3,n,av,bv,exx,ezz);
ElhL(:,l)=Elh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale_PSI;
computesISBdipoles;
computesIBoverlap;

EEc_hhL(:,:,l)  = EEc_hh;    % storing of the interband transition vs layer thickness
EEc_lhL(:,:,l)  = EEc_lh;    % storing of the interband transition vs layer thickness
EEhh_lhL(:,:,l) = EEhh_lh;   % storing of the interband transition vs layer thickness

EEc_cL(:,:,l)   = EEc_c;     % storing of the intersubband transition vs layer thickness
EEhh_hhL(:,:,l) = EEhh_hh;   % storing of the intersubband transition vs layer thickness
EElh_lhL(:,:,l) = EElh_lh;   % storing of the intersubband transition vs layer thickness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Display_IB_Results == 1
  PrintIBResults;
end
if Display_ISB_Results == 1
  PrintISBResults;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[10 50 1200 900]);
FS=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,1,'fontsize',FS)
hold on;grid on;

for i=1:length(EcL(:,1))
  plot(L,EcL(i,:),'o-')
  ll{i}=strcat('e',num2str(i));
end

legend(ll)
xlabel('Layer thickness (nm)')
ylabel('Ec (eV)')
title('electron in conduction band')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2,'fontsize',FS)
hold on;grid on;

for i=1:length(EhhL(:,1))
  plot(L,EhhL(i,:),'o-')
  ll{i}=strcat('hh',num2str(i));
end

legend(ll)
xlabel('Layer thickness (nm)')
ylabel('Ehh (eV)')
title('heavy hole in valence band')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,3,'fontsize',FS)
hold on;grid on;

for i=1:length(ElhL(:,1))
  plot(L,ElhL(i,:),'o-')
  ll{i}=strcat('lh',num2str(i));
end

legend(ll)
xlabel('Layer thickness (nm)')
ylabel('Elh (eV)')
title('light hole in valence band')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,4,'fontsize',FS)
hold on;grid on;

for i=1:length(EEc_hhL(1,:,1))
  plot(L,squeeze(EEc_hhL(i,i,:)),'o-')
  s1{i}=strcat('e',num2str(i),'-hh',num2str(i));
end
for i=1:length(EEc_lhL(1,:,1))
  plot(L,squeeze(EEc_lhL(i,i,:)),'o-')
  s2{i}=strcat('e',num2str(i),'-lh',num2str(i));
end

xlabel('Layer thickness (nm)')
ylabel('Ec-Eh (eV)')
title('Interband transition')
legend([s1 s2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n>1
  
subplot(2,3,5,'fontsize',FS)
hold on;grid on;
cc=1;
for i=1:length(EEc_cL(:,1,1))
  for j=1:length(EEc_cL(1,:,1))
    if j>i
      plot(L,squeeze(EEc_cL(i,j,:)),'o-')
      ss{cc}=strcat('e',num2str(i),'-e',num2str(j));
      cc=cc+1;
    end
  end
end

xlabel('Layer thickness (nm)')
ylabel('Ec-Ec (eV)')
title('Intersubband transition in CB')
legend(ss)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,6,'fontsize',FS)
hold on;grid on;
cc=1;
for i=1:length(EEhh_hhL(:,1,1))
  for j=1:length(EEhh_hhL(1,:,1))
    if j>i
      plot(L,squeeze(EEhh_hhL(i,j,:)),'o-')
      ss{cc}=strcat('hh',num2str(i),'-hh',num2str(j));
      plot(L,squeeze(EElh_lhL(i,j,:)),'o-')
      ss{cc+1}=strcat('lh',num2str(i),'-lh',num2str(j));
      cc=cc+2;
    end
  end
end

xlabel('Layer thickness (nm)')
ylabel('Eh-Eh (eV)')
title('Intersubband transition in VB')
legend(ss)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
