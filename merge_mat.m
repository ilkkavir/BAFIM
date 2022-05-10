function [success mergefile] = merge_mat(guisdapdir,append,delfiles)
%
% [success mergefile] = merge_mat(guisdapdir,append)
%
% Merge the guisdap matlab output files in a directory into one large matlab file.
% contents of each original file are stored to a separete structure array.
%
% INPUT:
%  guisdapdir  path to the guisdap output directory (result_path in guisdap)
%  append      if true, the data are appended to an existing file (if it exists) otherwise
%              a possible existing file is overwritten
%  delfiles    if true, the original matlab files are deleted
%
% OUTPUT:
%  success    0 if the merging was successful, otherwise see the error codes below
%  mergefile  full path to the merged output file
%
% Error codes
% 0 success
% 1 failed to list the guisdap files
% 2 failed to load a file
% 3 failed to save in the merged file
% 4 failed to delete an old file
%
% IV 2022
%

    success = 0;
    
    odir = strtrim(guisdapdir);
    nchar = length(odir);
    while odir(nchar)==filesep
        nchar = nchar-1;
    end
    odir = odir(1:nchar);
    if odir=='.'
        odir = pwd;
    end
    dirparts = strsplit(odir,filesep);
    
    mergefile = fullfile( odir , [char(dirparts(end)) '_merged.mat'] );

    % list the proper output files using guisdap getfilelist
    try
        guisdapfiles = getfilelist(fullfile(guisdapdir,filesep));
    catch
        success = 1;
        return;
    end

    for ifile = 1:length(guisdapfiles)
        try
            tmp = load(guisdapfiles(ifile).fname);
        catch
            success = 2;
            return
        end
        structname = sprintf('data%08d',guisdapfiles(ifile).file);
        eval([structname '=tmp;']);
        try
            if ifile==1
                if append
                    if exist(mergefile,'file')
                        save(mergefile,structname,'-append');
                    else
                        save(mergefile,structname);                    
                    end
                else
                    save(mergefile,structname);
                end
            else
                save(mergefile,structname,'-append');
            end
        catch
            success = 3;
            return
        end
        try
            if delfiles
                delete(guisdapfiles(ifile).fname);
            end
        catch
            success = 4;
        end
        disp(guisdapfiles(ifile).file)
    end
end

