function Pdat2mat(fnum, ifdel)
    % Transform .dat text file to .mat file.
    % Loading .dat file with matlab is time consuming.
	if nargin < 2
		ifdel = false;
	end
	if fnum(1) <= 0
		Fname = dir('particle*.dat');
		if ~isempty(Fname)
			for k = 1:numel(Fname)
				tic;
				P = loadp_f(Fname(k).name);
				filename = join([Fname(k).name(1:end-3), 'mat']);
				save(filename, 'P');
				fprintf([filename, ' Processed\n']);
				if ifdel
					delete(Fname(k).name);
				end
				toc;
			end
		end
	else
		for k = 1:numel(fnum)
			tic;
		    P = loadp_f(fnum(k));
		    filename = ['particle', num2str(fnum(k), '%04d')];
		    save(filename, 'P');
		    fprintf([filename, ' Processed\n']);
			if ifdel
				delete(join(['particle', num2str(fnum(k), '%04d'), '.dat']));
			end
			toc;
		end
	end
end