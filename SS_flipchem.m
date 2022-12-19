function SS = SS_flipchem(param,parinds,param_fit,covar,stdModErr,glat,glon,alt,fc)
%
% Sum-of-squares function to be minimized when adding flipchem information to the plasma parameter fits
%
% SS = SS_flipchem(param,parinds,param_fit,covar,stdModErr,date,glat,glon,alt)
%
% The flipchem python package must be installed (and the index files must be up to date)
%
% IV 2021
%


    % first copy the original fit results
    param_new = param_fit;

    % then replace the new values given in param
    param_new(parinds) = param;
    
    % convert the parameters to real units
    param_r = scaled_to_real(param_new);

    % call the model and update the composition
    ne = param_r(1);
    te = param_r(2)*param_r(3);
    ti = param_r(2);

    %    fc_comp = dirthe_flipchem(ne,te,ti,glat,glon,alt,fc);

    outputs = fc.get_point(glat,glon,alt,ne,te,ti);
    
    fc_comp = outputs{4}/ne;
    %    fc_comp = 1 - (outputs{5} + outputs{6})/ne;
    % % take the scaled parameter vector, we will change only the composition
    % param2 = param_new;
    % param2(6) = outputsm{4}/param_r(1);

    %    % select the parameters that were actually fitted
    %    idd = diag(covar)~=0;

    %    SS = (param_fit(idd)-param_new(idd))*inv(covar(idd,idd))*(param_fit(idd)-param_new(idd))' + (modErr-priorModErr)^2/stdModErr^2;
    SS = (param_fit(parinds)-param_new(parinds))*inv(covar(parinds,parinds))*(param_fit(parinds)-param_new(parinds))' + (param_new(6)-fc_comp)^2/stdModErr^2;

end