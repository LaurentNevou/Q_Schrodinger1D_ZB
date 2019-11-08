%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computes ISB dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Take care! Some people use meff inside the oscillator strenght f
% Actually, meff has sens in an infinite QW because there is a single mass value
% but NOT in multi-QW structure with various materials
% https://www.nextnano.com/nextnano3/tutorial/1Dtutorial_IntrabandTransitions.htm
% https://www.nextnano.com/nextnano3/tutorial/1Dtutorial_InGaAs_MQWs.htm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEc_c      = zeros(length(Ec),length(Ec));
z_dipole_c = zeros(length(Ec),length(Ec));
f_dipole_c = zeros(length(Ec),length(Ec));

for i=1:length(Ec)
  for j=1:length(Ec)
    if j>i
      EEc_c(i,j)      = Ec(j)-Ec(i);
      z_dipole_c(i,j) = abs(  trapz( z , psic(:,i).*z'.*psic(:,j) )  );
      f_dipole_c(i,j) = 2*m0/hbar^2 * ( Ec(j)-Ec(i) )* e * z_dipole_c(i,j)^2 ;
    end
  end
end

EEc_c      = EEc_c      + EEc_c.'     ;
z_dipole_c = z_dipole_c + z_dipole_c.';
f_dipole_c = f_dipole_c + f_dipole_c.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEhh_hh     = zeros(length(Ehh),length(Ehh));
z_dipole_hh = zeros(length(Ehh),length(Ehh));
f_dipole_hh = zeros(length(Ehh),length(Ehh));

for i=1:length(Ehh)
  for j=1:length(Ehh)
    if j>i
      EEhh_hh(i,j)     = abs(Ehh(j)-Ehh(i));
      z_dipole_hh(i,j) = abs(  trapz( z , psihh(:,i).*z'.*psihh(:,j) )  );
      f_dipole_hh(i,j) = 2*m0/hbar^2 * ( Ehh(j)-Ehh(i) )* e * z_dipole_hh(i,j)^2 ;
    end
  end
end

EEhh_hh      = EEhh_hh     + EEhh_hh.'    ;
z_dipole_hh  = z_dipole_hh + z_dipole_hh.';
f_dipole_hh  = f_dipole_hh + f_dipole_hh.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EElh_lh     = zeros(length(Elh),length(Elh));
z_dipole_lh = zeros(length(Elh),length(Elh));
f_dipole_lh = zeros(length(Elh),length(Elh));

for i=1:length(Elh)
  for j=1:length(Elh)
    if j>i
      EElh_lh (i,j)    = abs(Elh(j)-Elh(i));
      z_dipole_lh(i,j) = abs(  trapz( z , psilh(:,i).*z'.*psilh(:,j) )  );
      f_dipole_lh(i,j) = 2*m0/hbar^2 * ( Elh(j)-Elh(i) )* e * z_dipole_lh(i,j)^2 ;
    end
  end
end

EElh_lh      = EElh_lh     + EElh_lh.'    ;
z_dipole_lh  = z_dipole_lh + z_dipole_lh.';
f_dipole_lh  = f_dipole_lh + f_dipole_lh.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%