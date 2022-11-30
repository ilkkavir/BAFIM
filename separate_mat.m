function success = separate_mat(mergedfile,delfile)
%
% success = separate_mat(mergedfile,delfile)
%
% Separate data in a merged file into separate matlab files
%
% INPUT:
%  mergedfile  the merged output file or a directory containing one or more merged files
%  delfile     true to delete the merged file after separating the files
%
% OUTPUT:
%  success     0 if the process was successful, otherwise see error codes below
%
% Error codes
% 0 success
% 1 failed to read fieldnames from the input file
% 2 failed to load data from the input file
% 3 failed to save a mat file
% 4 failed to delete the input file
% 5 the input file is not a merged GUISDAP output file
%
% IV 2022
%
    success = 0;


    % if a directory was given as input
    if exist(mergedfile,'dir')
        % list matlab files in directory
        matfiles = dir(fullfile(mergedfile,'*.mat'));
        for k=1:length(matfiles)
            disp(matfiles(k).name)
            success(k) = separate_mat(fullfile(mergedfile,matfiles(k).name),delfile);
        end
        return
    end
    
    
    % the output directory
    dirparts = strsplit(mergedfile,filesep);
    if length(dirparts)==1
        odir = pwd;
    else
        odir = fileparts(mergedfile);
    end
    
    % field names in the merged file
    try
        fnames = sort(fieldnames(matfile(mergedfile)));
    catch
        success = 1;
        return
    end

    ismergedfile = false;
    for k=1:length(fnames)
        fname = char(fnames(k));
        if length(fname)==12
            if fname(1:4)=='data'
                % load contents of one separate mat file
                try
                    tmpdata = getfield(load(mergedfile,fname),fname);
                catch
                    success = 2;
                    return
                end
                outfile = fullfile(fileparts(mergedfile),[fname(5:12) '.mat']);
                try
                    save(outfile,'-struct','tmpdata');
                catch
                    success = 3;
                    return
                end
                ismergedfile = true;
                disp(outfile)
            end
        end
    end
    
    if delfile & ismergedfile
        try 
            delete(mergedfile)
        catch
            success = 4;
        end
    end

    if ~ismergedfile
        success = 5;
    end
    
end

