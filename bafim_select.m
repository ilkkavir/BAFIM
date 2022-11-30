function bafim_select( datadir , restype )
%
% bafim_select( datadir , restype )
%
% Switch between filtering and smoothing results in r_param and r_error variables in GUISDAP output files
%
% INPUT:
%   datadir  a directory containing bafim_smoother outputs
%
% UTPUT:
%   restype  'filter', 'smooth', 'rcorr', or 'r_cors'.
%            'rcorr' are the parameter profiles smoothed in range,
%            before adding the boundary conditions and process noise.
%
%
% IV 2020
%
    df = dir(fullfile(datadir,'*.mat'));

    nf  = length(df);

    for k=1:nf
        
        % read the data
        dfpath = fullfile(datadir,df(k).name);
        dd = load(dfpath);

        switch lower(restype(1:5))
          case 'filte'
            r_param = dd.r_param_filter;
            r_error = dd.r_error_filter;
          case 'smoot'
            r_param = dd.r_param_smooth;
            r_error = dd.r_error_smooth;
          case 'rcorr'
            r_param = dd.r_param_rcorr;
            r_error = dd.r_error_rcorr;
          case 'rcors'
            r_param = dd.r_param_rcorr_smooth;
            r_error = dd.r_error_rcorr_smooth;
        otherwise
            error("unknown restype, must be either 'filter' or 'smooth'")
        end

        save(dfpath,'r_param','r_error','-append');

        fprintf("\r %s",dfpath)
        
    end

    
end

