%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Print IB Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('')
display('===================================================')
display('Interband Results:')
display('===================================================')
display('')

for i=1:length(Ec)
  for j=1:length(Ehh)
    
    display(strcat(...
    'e',num2str(i),'-hh',num2str(j),' = ',num2str( EEc_hh(i,j) ,'%.3f' ),'eV = ',num2str(h*c/e/EEc_hh(i,j)*1e9,'%.0f'),'nm = '...
    ,num2str(EEc_hh(i,j)/h*e*1e-12,'%.0f'),'THz = ',num2str(EEc_hh(i,j)*e/h/c*1e-2,'%.0f'),'cm-1; |<e',num2str(i),'|hh',num2str(j),'>|='...
    ,num2str(IBoverlap_ehh(i,j),'%.3f'),'; |<e',num2str(i),'|hh',num2str(j),'>|2=',num2str(IBoverlap_ehh(i,j)^2,'%.3f')...
    ) )
    
  end
end

display('')

for i=1:length(Ec)
  for j=1:length(Elh)
    
    display(strcat(...
    'e',num2str(i),'-lh',num2str(j),' = ',num2str( EEc_lh(i,j) ,'%.3f' ),'eV = ',num2str(h*c/e/EEc_lh(i,j)*1e9,'%.0f'),'nm = '...
    ,num2str(EEc_lh(i,j)/h*e*1e-12,'%.0f'),'THz = ',num2str(EEc_lh(i,j)*e/h/c*1e-2,'%.0f'),'cm-1; |<e',num2str(i),'|lh',num2str(j),'>|='...
    ,num2str(IBoverlap_elh(i,j),'%.3f'),'; |<e',num2str(i),'|lh',num2str(j),'>|2=',num2str(IBoverlap_elh(i,j)^2,'%.3f')...
    ) )
    
  end
end

display('')

for i=1:length(Ehh)
  for j=1:length(Elh)
    
    display(strcat(...
    'hh',num2str(i),'-lh',num2str(j),' = ',num2str( EEhh_lh(i,j) ,'%.3f' ),'eV = ',num2str(h*c/e/EEhh_lh(i,j)*1e6,'%.1f'),'um = '...
    ,num2str(EEhh_lh(i,j)/h*e*1e-12,'%.1f'),'THz = ',num2str(EEhh_lh(i,j)*e/h/c*1e-2,'%.0f'),'cm-1; |<hh',num2str(i),'|lh',num2str(j),'>|='...
    ,num2str(IBoverlap_hhlh(i,j),'%.3f'),'; |<hh',num2str(i),'|lh',num2str(j),'>|2=',num2str(IBoverlap_hhlh(i,j)^2,'%.3f')...
    ) )
    
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
