function bafim_smoother( datadir , mergedfile )
%
% bafim_smoother( datadir )
%
% An RTS-smoother for GUISDAP BAFIM output.
%
%
% INPUT:
%   datadir     a GUISDAP output directory with BAFIM fit results
%   mergedfile  true if the input data merged guisdap output files, false for normal guisdap files
%
% OUTPUT:
%             none, the smoothed results are appended to the original output files
%
%
%
% IV 2020
%
    global path_GUP

    persistent r_param_smooth_next r_error_smooth_next r_apriori_next r_apriorierror_next

    if nargin < 2
        mergedfile = false;
    end
    
    % if isdir(datadir)
    %     df = dir(fullfile(datadir,'*.mat'));
    %     nf  = length(df);
    %     merfedfile = false;
    % else
    %     fnames = sort(fieldnames(matfile(datadir)));
    %     l = 1;
    %     flist = [];
    %     for k=1:length(fnames)
    %         fname = char(fnames(k));
    %         if length(fname)==12
    %             if fname(1:4)=='data'
    %                 flist(l).file = str2num(fname(5:12));
    %                 flist(l).fname = fullfile(fileparts(datadir),[fname(5:12) '.mat']);
    %                 l = l+1;
    %             end
    %         end
    %     end
    %     nf = length(flist); 
    %     mergedfile = true;
    % end

    if mergedfile
        % file nams of the merged files must start with "GUISDAP"
        df = dir(fullfile(datadir,'GUISDAP*.mat'));
        l = 1;
        flist = [];
        for ifile=1:length(df)
            % the file names should be properly ordered by dir
            fnames = sort(fieldnames(matfile(fullfile(df(ifile).folder,df(ifile).name))));
            for k=1:length(fnames)
                fname = char(fnames(k));
                if length(fname)==12
                    if fname(1:4)=='data'
                        flist(l).file = str2num(fname(5:12));
                        flist(l).fname = fullfile(datadir,[fname(5:12) '.mat']);
                        flist(l).mergedfname = df(ifile).name;
                        l = l+1;
                    end
                end
            end
        end
        nf = l-1;
    else
        df = dir(fullfile(datadir,'*.mat'));
        nf  = length(df);
    end

    % physically reasonable limits for the plasma parameters
    paramlims = [1e6 1 .01 1 -2e4 -.01 ; 1e14 2e4 100 1e9 2e4 1.01];

    for k=nf:-1:1
        if k==nf
            addpath(fullfile(path_GUP,'init'))
        end

        % read the data
        if ~mergedfile
            dfpath = fullfile(datadir,df(k).name);
            dd = load(dfpath);
        else
            fname = ['data' flist(k).fname((end-11):(end-4))];
            dd = getfield(load(fullfile(datadir,flist(k).mergedfname),fname),fname);
        end

        % number of height gates (actually, the present version cannot handle changes in nhei)
        nhei = length(dd.r_h);

        % The unsmoothed data will be written in r_param_filter and r_error_filter, and 
        % the smoothed ones in both r_param_smooth & r_error_smooth, and r_param & r_error.
        % Use r_param_filter & r_error_filter if they exists, otherwise read r_param & r_error.
        % Then we do not have problems if the smoother is accidentially run more than once


        if isfield(dd,'r_param_filter')
            r_param_filter = dd.r_param_filter;
        else
            r_param_filter = dd.r_param;
        end
        if isfield(dd,'r_error_filter')
            r_error_filter = dd.r_error_filter;
        else
            r_error_filter = dd.r_error;
        end

        if k < nf

            r_param_filter_s = real_to_scaled( r_param_filter );
            r_error_filter_s = real_to_scaled( r_error_filter);
            r_param_smooth_next_s = real_to_scaled( r_param_smooth_next);
            r_error_smooth_next_s = real_to_scaled( r_error_smooth_next);
            r_apriori_next_s = real_to_scaled( r_apriori_next );
            r_apriorierror_next_s = real_to_scaled( r_apriorierror_next );

            P_k = zeros(nhei*5,nhei*5);
            P_k1_smooth = zeros(nhei*5,nhei*5);
            m_k = NaN(nhei*5,1);
            m_k1_smooth = NaN(nhei*5,1);
            m_k1_pred = NaN(nhei*5,1);
            P_k1_pred = zeros(nhei*5,nhei*5);
            for ihei = 1:nhei
                fitcov = vec2covm(r_error_filter_s(ihei,:));
                P_k( ((0:4)*nhei+ihei) , ((0:4)*nhei+ihei) ) = fitcov([1 2 ...
                                    3 5 6],[1 2 3 5 6]);
                m_k( ((0:4)*nhei)+ihei ) = r_param_filter_s(ihei,[1 2 3 5 6]);

                fitcov = vec2covm(r_error_smooth_next_s(ihei,:));
                P_k1_smooth( ((0:4)*nhei+ihei) , ((0:4)*nhei+ihei) ) = fitcov([1 2 ...
                                    3 5 6],[1 2 3 5 6]);
                m_k1_smooth( ((0:4)*nhei)+ihei ) = r_param_smooth_next_s(ihei,[1 2 3 5 6]);

                m_k1_pred( ((0:4)*nhei)+ihei ) = r_apriori_next_s(ihei,[1 2 3 5 6]);

                P_k1_pred( ((0:4)*nhei+ihei) , ((0:4)*nhei+ihei) ) = diag(r_apriorierror_next_s(ihei,[1 2 3 5 6]).^2);
            end

            m_k_smooth = m_k + dd.BAFIM_G * ( m_k1_smooth - m_k1_pred);
            P_k_smooth = P_k + dd.BAFIM_G * ( P_k1_smooth - P_k1_pred ) * dd.BAFIM_G';

            r_param_smooth_s = real_to_scaled(r_param);
            r_param_smooth_s(:,1) = m_k_smooth(1:nhei);
            r_param_smooth_s(:,2) = m_k_smooth((nhei+1):(2*nhei));
            r_param_smooth_s(:,3) = m_k_smooth((2*nhei+1):(3*nhei));
            r_param_smooth_s(:,5) = m_k_smooth((3*nhei+1):(4*nhei));
            r_param_smooth_s(:,6) = m_k_smooth((4*nhei+1):(5*nhei));

            r_error_smooth_s = real_to_scaled(r_error);
            var_k_smooth = diag(P_k_smooth);
            r_error_smooth_s(:,1) = sqrt(var_k_smooth(1:nhei));
            r_error_smooth_s(:,2) = sqrt(var_k_smooth((nhei+1):(2*nhei)));
            r_error_smooth_s(:,3) = sqrt(var_k_smooth((2*nhei+1):(3*nhei)));
            r_error_smooth_s(:,5) = sqrt(var_k_smooth((3*nhei+1):(4*nhei)));
            r_error_smooth_s(:,6) = sqrt(var_k_smooth((4*nhei+1):(5*nhei)));
            

            r_param_smooth = scaled_to_real(r_param_smooth_s);
            r_error_smooth = scaled_to_real(r_error_smooth_s);
        else
            r_param_smooth = r_param_filter;
            r_error_smooth = r_error_filter;
        end

        % make sure that the smoothed values are within paramlims, this should rarely have any effect
        for ihei=1:nhei
            r_param_smooth(ihei,1:6) = max(r_param_smooth(ihei,1:6),paramlims(1,:));
            r_param_smooth(ihei,1:6) = min(r_param_smooth(ihei,1:6),paramlims(2,:));
        end

        % points with std <= 1e-3 are points fixed to model values, put the model values in these points
        r_param_smooth(r_error_filter(:,1:6)<=1e-3) = r_param_filter(r_error_filter(:,1:6)<=1e-3);
        
        r_param_smooth_next = r_param_smooth;
        r_error_smooth_next = r_error_smooth;
        r_apriori_next = dd.r_apriori;
        r_apriorierror_next = dd.r_apriorierror;


        r_param = r_param_smooth;
        r_error = r_error_smooth;

        if ~mergedfile
            save(dfpath,'r_param','r_param_smooth','r_param_filter','r_error','r_error_smooth','r_error_filter','-append');
            fprintf("\r %s",dfpath)
        else

            % this is slow even with one-hour files. Should we first save individual files and then merge them?
            % -- may work otherwise, but would need to make sure that too many files do not exist at the same time
            dd.r_param = r_param;
            dd.r_param_smooth = r_param_smooth;
            dd.r_param_filter = r_param_filter;
            dd.r_error = r_error;
            dd.r_error_smooth = r_error_smooth;
            dd.r_error_filter = r_error_filter;

            structname = sprintf('data%08d',flist(k).file);
            eval([structname '=dd;']);
            save(fullfile(datadir,flist(k).mergedfname),structname,'-append')

            % dfpath = fullfile(datadir,flist(k).fname);
            % save(dfpath,'-struct','dd')

            
            fprintf("\r %08d",flist(k).file)
            %            fprintf("\r %s",dfpath)

        end

        
        
    end


    % if mergedfile
    %     merge_mat(datadir,false,true)
    % end

end