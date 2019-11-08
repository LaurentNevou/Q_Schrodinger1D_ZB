%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Computes IB wavefunctions overlap %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% https://www.nextnano.de/nextnano3/tutorial/1Dtutorial_OpticalTransitions.htm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IBoverlap_ehh  = zeros(length(Ec) ,length(Ehh));
IBoverlap_elh  = zeros(length(Ec) ,length(Elh));
IBoverlap_hhlh = zeros(length(Ehh),length(Elh));


for i=1:length(Ec)
  for j=1:length(Ehh)
    
    EEc_hh(i,j) = Ec(i)-Ehh(j);
    IBoverlap_ehh(i,j) = abs(  trapz( z , psic(:,i).*psihh(:,j) )  );
    
  end
end

for i=1:length(Ec)
  for j=1:length(Elh)
    
    EEc_lh(i,j) = Ec(i)-Elh(j);
    IBoverlap_elh(i,j) = abs(  trapz( z , psic(:,i).*psilh(:,j) )  );
    
  end
end

for i=1:length(Ehh)
  for j=1:length(Elh)
    
    EEhh_lh(i,j) = abs(Ehh(i)-Elh(j));
    IBoverlap_hhlh(i,j) = abs(  trapz( z , psihh(:,i).*psilh(:,j) )  );
    
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%