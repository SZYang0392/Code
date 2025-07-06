function Fdat2mat(fnum)
    % Transform .dat text file to .mat file.
    % Loading .dat file with matlab is time consuming.
    
    for k = 1:numel(fnum)
        F = loadfield(fnum(k));
        filename = ['field', num2str(fnum(k), '%05d')];
        save(filename, 'F');
        fprintf([filename, ' Processed\n']);
    end
end