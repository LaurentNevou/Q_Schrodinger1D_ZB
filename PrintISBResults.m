%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Print ISB Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('')
display('===================================================')
display('Intersubband Results:')
display('===================================================')
display('')

for i=1:length(Ec)
  for j=1:length(Ec)
    if j>i
        display(strcat(...
        'e',num2str(i),'-e',num2str(j),' = ',num2str( EEc_c(i,j),'%.3f' ),'eV;   z'...
        ,num2str(i),'-',num2str(j),' = ',num2str( z_dipole_c(i,j)*1e9,'%.3f' ),'nm;   f'...
        ,num2str(i),'-',num2str(j),' = ',num2str( f_dipole_c(i,j),'%.3f' ) ...
        )  )
    end
  end
end

display('')

for i=1:length(Ehh)
  for j=1:length(Ehh)
    if j>i
        display(strcat(...
        'hh',num2str(i),'-hh',num2str(j),' = ',num2str( EEhh_hh(i,j),'%.3f' ),'eV;   z'...
        ,num2str(i),'-',num2str(j),' = ',num2str( z_dipole_hh(i,j)*1e9,'%.3f' ),'nm;   f'...
        ,num2str(i),'-',num2str(j),' = ',num2str( f_dipole_hh(i,j),'%.3f' ) ...
        )  )
    end
  end
end

display('')

for i=1:length(Elh)
  for j=1:length(Elh)
    if j>i
        display(strcat(...
        'lh',num2str(i),'-lh',num2str(j),' = ',num2str( EElh_lh(i,j),'%.3f' ),'eV;   z'...
        ,num2str(i),'-',num2str(j),' = ',num2str( z_dipole_lh(i,j)*1e9,'%.3f' ),'nm;   f'...
        ,num2str(i),'-',num2str(j),' = ',num2str( f_dipole_lh(i,j),'%.3f' ) ...
        )  )
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%