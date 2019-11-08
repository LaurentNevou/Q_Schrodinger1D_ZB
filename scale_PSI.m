%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Scaling and shifting the wavefunctions %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSIc=[];PSIhh=[];PSIlh=[];

for i=1:length(Ec)
    PSIc(:,i)=abs(psic(:,i)).^2/max(abs(psic(:,i)).^2)*ScF + Ec(i); % normalisation for the plotting
end
for i=1:length(Ehh)
    PSIhh(:,i)=abs(psihh(:,i)).^2/max(abs(psihh(:,i)).^2)*ScF + Ehh(i); % normalisation for the plotting
end
for i=1:length(Elh)
    PSIlh(:,i)=abs(psilh(:,i)).^2/max(abs(psilh(:,i)).^2)*ScF + Elh(i); % normalisation for the plotting
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%