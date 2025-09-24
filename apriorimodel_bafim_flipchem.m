function [apriori2,apriorierror2] = apriorimodel_bafim_flipchem(apriori,apriorierror,heights,fit_altitude)
%
%
% [apriori2,apriorierror2] = apriorimodel_bafim_flipchem(apriori,apriorierror,heights)
%
% Bayesian filtering in time and correlation priors in range direction for Ne, Ti, Tr, Vi, and composition.
% Ion composition is updated with the flipchem model in each prediction step.
%
% NOTICE: the whole BAFIM package, including this function, ionomodel_bafim, and the modified ionomodel-function must be in the MATLAB search path
%
%
%         Some GUISDAP input paramters will need to be correctly set to use BAFIM (you can write these in the "extra" box of the "guisdap for dummies" window):
%         iono_model='bafim'
%         a_phasepush=0
%         a_satch.do=0
%         fit_altitude(1:6,1:4) = [h1Ne h2Ne hsNe tsNe ; h1Ti h2Ti hsTi tsTi ; h1Tr h2Tr hsTr tsTr ; 0 0 0 0 ; h1Vi h2Vi hsVi tsVi ; h1O+ h2O+ hsO+ tsO+ ],
%         where
%            h1Ne is the lowest altitude where Ne is fitted
%            h2Ne is the highest altitude  where Ne is fitted
%            hsNe is additional scaling for Ne correlation length (see below)
%            tsNe is scale of the process noise standard deviation for Ne (see below)
%            h1Ti, h2Ti, etc. are the corresponding limits and scales for the other parameters
%
%
%            hsNe, hsTi, hsTr, hsVi and hsO+ are additional scales for correlation length.
%            The final correlation is the product of the additional scaling factor and plasma scale height (calculated from IRI parameters)
%
%            tsNe, tsTi, tsTr, tsVi and tsO+ are scales of the process noise standard deviation
%            The final process noise standard devition is the product of the scale factor and square root of the time step (in seconds)
%
%        A combination of limits and scaling factors that works reasonably well with the default GUISDAP altitude grid and ~5s time resolution is:
%        fit_altitude(1:6,1:4) = [0 Inf 0.1 2e11 ; 80 Inf 0.3 10 ; 103 Inf 0.3 0.05 ; 0 0 0 0 ; 80 Inf 0.2 2.5 ; 130 350 .2 .003]
%
%        Notice that the system is grid dependent. In general, larger process noise variances and correlation lengths are needed for higher range resolution.
%        The optimal tuning depends on SNR in a non-trivial way. Grid-independnce and automatic tuning with SNR might be implemented in future...
%
%        The option a_phasepush=0 disables phase pushing correction in GUISDAP, and a_stach.do disables the satellite echo detection. The phase pushing
%        correction and the satelite echo removal were unfortunately implemented in a way that do not fit together fith this version of bafim.
%
%
%
% INPUT:
%        apriori       initial prior parameters as returned ionomodel_bafim
%        apriorierror  prior errors calculated in ionoomdel. The matrix is not used in this version of bafim.
%        heights       centre points of height gates
%        fit_altitude  the fit_altitude matrix used internally in GUISDAP
%
% OUTPUT:
%        apriori2      predicted parameters of the Bayesian filter. These are used as a prior in the iterative fit.
%        apriorierror2 standard deviations of the predicted parameters.
% 
%
% See also: bafim_flipchem_smoother
%
% IV 2018-2025
%
    
    
% Contents of some global GUISDAP variables for reference
%
% r_param   plasma parameter fit results from the previous time
%           step
% r_error   error estiamtes for the plasma parameters
% r_res     residuals of the fits in previous time step
% r_status  fit status for each height gate in the previous time
%           step (0 if the fit is ok)
% d_time    timestap from the raw data
    
    global r_param r_error r_res r_status d_time path_GUP result_path v_Boltzmann v_amu p_XMITloc
    global v_lightspeed v_Boltzmann v_epsilon0 v_elemcharge sc_angle k_radar

    % merge outputs in one large file if true, set to false for the normal guisdap output files
    merge_output_files = false;

    % flipchem model error standar deviation
    flipchem_modErr_std = .1;

    % step length in finite difference gradient calculations
    dx = .001;
    
    % We need to pass some variables from one timestep to another
    persistent d_time_prev aprioriprev apriorierrorprev filename istep tfromstart

    % make a  copy of the original GUISDAP fit results
    r_param_orig = r_param;
    r_error_orig = r_error;

    % Lowest & highest altitudes for fits from the input array
    hlimNe = fit_altitude(1,1:2);
    hlimTi = fit_altitude(2,1:2);
    hlimTr = fit_altitude(3,1:2);
    hlimVi = fit_altitude(5,1:2);
    hlimOp = fit_altitude(6,1:2);

    % number of height gates
    nhei = length(heights);
    
    % the lowest gate must be from a model, except in remote analysis
    if nhei > 3
        hlimTi(1) = max(hlimTi(1),(min(heights)+.01));
        hlimTr(1) = max(hlimTr(1),(min(heights)+.01));
        hlimVi(1) = max(hlimVi(1),(min(heights)+.01));
        hlimOp(1) = max(hlimOp(1),(min(heights)+.01));
        % make sure that we have a model composition at the upper boundary
        hlimOp(2) = min(hlimOp(2),(max(heights)-.01));
    end

    % read the correlation length and process noise scales from the input array
    hsNe = fit_altitude(1,3);
    hsTi = fit_altitude(2,3);
    hsTr = fit_altitude(3,3);
    hsVi = fit_altitude(5,3);
    hsOp = fit_altitude(6,3);

    tsNe = fit_altitude(1,4);
    tsTi = fit_altitude(2,4);
    tsTr = fit_altitude(3,4);
    tsVi = fit_altitude(5,4);
    tsOp = fit_altitude(6,4);

    % limits for plasma parameters, these are more strict than those in fit_altitude(:,7:8)
    paramlims = [1e9 50 .1 1 -1e4 0 ; 1e13 1e4 10 1e9 1e4 1];
    
    % r_param is empty on the first iteration. Initialize d_time_prev, add the guisdap init directory to the path, and copy this file to the result directory.
    if size(r_param,1)==0
        d_time_prev = d_time(1,:);
        addpath([path_GUP 'init']);
        istep = 1;
        tfromstart = 0;
    end
    if istep==2
        % copy the prior model in the output directory to help in development work
        % copy only on the second step to be sure that result_path has been properly parsed
        copyfile(which('apriorimodel_bafim_flipchem.m'),fullfile(result_path,'apriorimodel_bafim_flipchem.m'));
    end

    
    % Time step length
    tstep = seconds(datetime(d_time(2,:))-datetime(d_time_prev));
    d_time_prev = d_time(2,:);
    
    % The final correlation lengths are products hsXX*H, where H is the plasma scale height as calculated from IRI parameters.
    H = v_Boltzmann .* apriori(:,2).*(1+apriori(:,3))./2 ./ ( v_amu .* ( apriori(:,6)*16 + 30.5*(1-apriori(:,6)) ) .* 9.82.*(6372./(6372+heights)).^2 );
    hsAlt = H/1000;

    % The correlation lengths must scale with time step duration to make the system (somewhat) independent of the step length
    hsAlt = hsAlt * sqrt(tstep);
    
    if nhei > 1
        % Approximate widths of the height gates
        dheights = diff(heights);
        dheights = [dheights(1) ; dheights];
    end
        
    % Copy the default prior model
    apriori2 = apriori;
    apriorierror2 = apriorierror;

    % an array for k * Debye length
    kdeb = [];

    % counters to re-initialize the filter after data gaps longer than 600 s
    if tstep > 600
        tfromstart = 0;
    else
        tfromstart = tfromstart + tstep;
    end

    % take prior from IRI if exp start was less than 2 minutes ago to avoid problems with radar system instabilities
    if tfromstart > 120

    
        % Replace unreasonable values with the previous predictions. 
        for hind = 1:nhei
            
            okfit = true;

            % obviously erroneous fits
            if ~any(r_status(hind)==[0 3]) | any(isnan(r_param(hind,1:6))) | any(isnan(r_error(hind,:))) |...
                    r_res(hind,1)>100 | any(r_param(hind,1:5)<paramlims(1,1:5)) | ...
                    any(r_param(hind,1:5)>paramlims(2,1:5))
                okfit = false;
            end

            % try to exclude large jumps due to weak space debris echoes etc.
            SSlim = 25;
            if heights(hind)>hlimNe(1) & heights(hind)<hlimNe(2)
                if ( (r_param(hind,1)-aprioriprev(hind,1))^2 / apriorierrorprev(hind,1)^2  ) > SSlim
                    okfit = false;
                end
            end
            if heights(hind)>hlimTi(1) & heights(hind)<hlimTi(2)
                if ( (r_param(hind,2)-aprioriprev(hind,2))^2 / apriorierrorprev(hind,2)^2  ) > SSlim
                    okfit = false;
                end
                if r_param(hind,2) < .3*apriori(hind,2)
                    okfit = false;
                end
            end
            if heights(hind)>hlimTr(1) & heights(hind)<hlimTr(2)
                if ( (r_param(hind,3)-aprioriprev(hind,3))^2 / apriorierrorprev(hind,3)^2  ) > SSlim
                    okfit = false;
                end
            end
            if heights(hind)>hlimVi(1) & heights(hind)<hlimVi(2)
                if ( (r_param(hind,5)-aprioriprev(hind,5))^2 / apriorierrorprev(hind,5)^2  ) > SSlim
                    okfit = false;
                end
            end
            if heights(hind)>hlimOp(1) & heights(hind)<hlimOp(2)
                if ( (r_param(hind,6)-aprioriprev(hind,6))^2 / apriorierrorprev(hind,6)^2  ) > SSlim
                    okfit = false;
                end
            end
            if ~okfit
                r_param(hind,1:6) = aprioriprev(hind,1:6);
                r_error(hind,1:6) = apriorierrorprev(hind,1:6);
                r_error(hind,7:end) = 0;
            end

            % set r_status to zero if it was 3 only because O+ fraction is slightly above 1
            if r_status(hind)==3 & r_param(hind,6)>1 & r_param(hind,6)<1.1
                r_status(hind) = 0;
            end
            % set r_status to zero if it was 3 only because O+ fraction is slightly below 0
            if r_status(hind)==3 & r_param(hind,6)<0 & r_param(hind,6)>-.1
                r_status(hind) = 0;
            end
        end

        % exclude points that suffer from finite Debye length on the topside
        if nhei > 5
            iweak = nhei+1;

            % Te and Ti, take max of model and fit
            tetmp = max(r_param(:,2).*r_param(:,3),apriori(:,2).*apriori(:,3));
            titmp = max(r_param(:,2),apriori(:,2));

            % Debye length
            debye = sqrt((v_epsilon0*v_Boltzmann/v_elemcharge^2)./(r_param(:,1).*(1./tetmp + 1./titmp)));

            % remove points where the debye length is a signficant fraction of the radar wavelength
            kdeb = debye .* k_radar(1);
            while any( kdeb((iweak-5):(iweak-1)) > .6)
                iweak = iweak - 1;
                if iweak == 6
                    break
                end
            end
            if iweak <= nhei
                for hind = iweak:nhei
                    % reject the point if Te/Ti fit was attempted, but keep in the bottom part
                    if heights(hind)>hlimTr(1) & heights(hind)<hlimTr(2)
                        r_param(hind,1:6) = aprioriprev(hind,1:6);
                        r_error(hind,1:6) = apriorierrorprev(hind,1:6);
                        r_error(hind,7:end) = 0;
                        r_status(hind) = 6;
                    end
                end
            end
            %            disp([iweak heights(iweak-1)])
        end

        % make a copy of the cleaned parameters
        r_param_orig_cleaned = r_param;
        r_error_orig_cleaned = r_error;
        % in scaled units
        r_param_orig_cleaned_s = real_to_scaled(r_param);
        r_error_orig_cleaned_s = real_to_scaled(r_error);
        
        % Convert to scaled units
        r_param_s = real_to_scaled(r_param);
        r_error_s = real_to_scaled(r_error);
        apriori_s = real_to_scaled(apriori);
        apriorierror_s = real_to_scaled(apriorierror);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update composition with flipchem %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fit all parameters with an additional penalty for deviations from the flipchem composition
        %
        % call flipchem directly from the python module
        %


        % initialize flipchem for the measurement time
        timetmp = datetime(d_time(2,:));
        pydate = py.datetime.datetime(int32(timetmp.Year),int32(timetmp.Month),int32(timetmp.Day),int32(timetmp.Hour),int32(timetmp.Minute));
        glat = p_XMITloc(1);
        glon = p_XMITloc(2);
        fc = py.flipchem.Flipchem(pydate);
        O2p = r_param(:,1).*NaN;
        NOp = r_param(:,1).*NaN;
        NO  = r_param(:,1).*NaN;

        % fms options...
        fms_opts = optimset('fminsearch');
        % fms_opts.Display = 'off';
        % fms_opts.MaxFunEvals=1e4;
        % fms_opts.MaxIter=1e6;
        % fms_opts.TolFun=1e-8;
        % fms_opts.TolX=1e-8;

        % update the parameters using flipchem at each altitude
        for hind = 1:nhei

            % altitude
            galt = heights(hind);

            % correct only if we fit composition in this gate
            if heights(hind)>hlimOp(1) & heights(hind)<hlimOp(2)

                r_param_s_copy = r_param_s(hind,:);
                r_param_copy = r_param(hind,:);
                r_error_s_copy = r_error_s(hind,:);
                r_error_copy = r_error(hind,:);
                % skip the composition update if something fails
                try
                    
                    % the error covariance matrix from previous time step (in scaled units)
                    covmat = vec2covm(r_error_s(hind,:));
                    
                    % the parameters that were actually fitted
                    fitpar = diag(covmat) ~= 0;
                    ifitpar = find(fitpar);

                    % the model error with additional scaling by height above 150 km (the std is doubled at every 350 km, because the model is less reliable at heigher altitudes)
                    fcStd = max(1,exp((heights(hind)-150)/(350/log(2)))) * flipchem_modErr_std;

                    % a new fit with penalty from deviation from flipchem
                    [x,fval,exitflag,output] = fminsearch( @(p) SS_flipchem( p , fitpar , r_param_orig_cleaned_s(hind,:) , covmat, fcStd , glat , glon , galt , fc ) , r_param_orig_cleaned_s(hind,fitpar) , fms_opts);
                    
                    % update the parameter vector
                    r_param_s(hind,fitpar) = x;
                    r_param(hind,:) = scaled_to_real(r_param_s(hind,:));
                    
                    
                    %%%%%%%%%%%%%%%%%%%%
                    % error estimation %
                    %%%%%%%%%%%%%%%%%%%%
                    
                    % number of parameters we actually fit
                    nparam = sum(fitpar);
                    
                    % initialize the Hessian matrix
                    Hess = zeros(nparam,nparam);
                    
                    % sum of squares at the iteration convergence point
                    SS0 = SS_flipchem( r_param_s(hind,fitpar) , fitpar , r_param_orig_cleaned_s(hind,:) , covmat, fcStd , glat , glon , galt , fc);
                    
                    % Diagonal of the hessian
                    for k=1:nparam
                        
                        % +dx
                        r_param_tmp = r_param_s(hind,:);
                        r_param_tmp(ifitpar(k)) = r_param_tmp(ifitpar(k)) + dx;
                        Hess(k,k) = Hess(k,k) + SS_flipchem( r_param_tmp(fitpar) , fitpar , r_param_orig_cleaned_s(hind,:) , covmat, fcStd , glat , glon , galt , fc); 
                        
                        % -dx
                        r_param_tmp = r_param_s(hind,:);
                        r_param_tmp(ifitpar(k)) = r_param_tmp(ifitpar(k)) - dx;
                        Hess(k,k) = Hess(k,k) + SS_flipchem( r_param_tmp(fitpar) , fitpar , r_param_orig_cleaned_s(hind,:) , covmat, fcStd , glat , glon , galt , fc);
                        
                        % the centre and scaling
                        Hess(k,k) = (Hess(k,k) - 2*SS0) / dx^2;
                        
                    end
                    
                    % the off-diagonal elements
                    for k=1:(nparam-1)
                        
                        for l=(k+1):nparam
                            
                            % +dx+dx
                            r_param_tmp = r_param_s(hind,:);
                            r_param_tmp(ifitpar(k)) = r_param_tmp(ifitpar(k)) + dx;
                            r_param_tmp(ifitpar(l)) = r_param_tmp(ifitpar(l)) + dx;
                            Hess(k,l) = Hess(k,l) + SS_flipchem( r_param_tmp(fitpar) , fitpar , r_param_orig_cleaned_s(hind,:) , covmat, fcStd , glat , glon , galt , fc); 
                            
                            % +dx-dx
                            r_param_tmp = r_param_s(hind,:);
                            r_param_tmp(ifitpar(k)) = r_param_tmp(ifitpar(k)) + dx;
                            r_param_tmp(ifitpar(l)) = r_param_tmp(ifitpar(l)) - dx;
                            Hess(k,l) = Hess(k,l) - SS_flipchem( r_param_tmp(fitpar) , fitpar , r_param_orig_cleaned_s(hind,:) , covmat, fcStd , glat , glon , galt , fc); 
                            
                            % -dx+dx
                            r_param_tmp = r_param_s(hind,:);
                            r_param_tmp(ifitpar(k)) = r_param_tmp(ifitpar(k)) - dx;
                            r_param_tmp(ifitpar(l)) = r_param_tmp(ifitpar(l)) + dx;
                            Hess(k,l) = Hess(k,l) - SS_flipchem( r_param_tmp(fitpar) , fitpar , r_param_orig_cleaned_s(hind,:) , covmat, fcStd , glat , glon , galt , fc); 
                            
                            % -dx-dx
                            r_param_tmp = r_param_s(hind,:);
                            r_param_tmp(ifitpar(k)) = r_param_tmp(ifitpar(k)) - dx;
                            r_param_tmp(ifitpar(l)) = r_param_tmp(ifitpar(l)) - dx;
                            Hess(k,l) = Hess(k,l) + SS_flipchem( r_param_tmp(fitpar) , fitpar , r_param_orig_cleaned_s(hind,:) , covmat, fcStd , glat , glon , galt , fc); 
                            
                            Hess(k,l) = Hess(k,l)/(4*dx^2);
                            Hess(l,k) = Hess(k,l);
                            
                        end
                    end
                    
                    % regularize a bit if there are unrealistic values (this happens very rarely but would crash the analysis)
                    while any(diag(Hess)<=0)
                        %                    disp('negative in Hess matrix diagonal')
                        Hess = Hess + eye(nparam);
                    end
                    
                    % inverse of the Hessian matrix
                    invhess = inv(Hess);
                    
                    % The diagonal should be positive, regularize if it is not (again, very rare)
                    while(any(diag(invhess)<=0))
                        %                    disp('negative in covariance matrix diagonal')
                        Hess = Hess + eye(nparam);
                        invhess = inv(Hess);
                    end
                    
                    % covariance matrix of all parameters
                    covmat2 = covmat;
                    
                    % error covariance = inv(2*Hess), the remaining elements of covmat2 are not affected
                    covmat2(fitpar,fitpar) = invhess*2;
                    
                    % convert the covariance matrix into the vector format used in guisdap
                    r_error_s(hind,:) = covm2vec(covmat2);

                    % the vector format in physical units
                    r_error(hind,:) = scaled_to_real(r_error_s(hind,:));

                    % final check that everything is ok after the flipchem fit
                    if any(any(isnan(covmat2)))
                        error('NaN in covariance matrix, skipping the flipchem fit.');
                    end
                    if any(any(diag(covmat2)<0))
                        error('Negative variance, skipping the flipchem fit.');
                    end
                    if any(isnan(r_error(hind,:)))
                        error('NaN standard deviation, skipping the flipchem fit.');
                    end
                    if any(r_error(hind,1:6)<0)
                        error('Negative standard deviation, skipping the flipchem fit.');
                    end

                    % if the flipchem fit fails we continue with the guisdap fit results
                catch me
                    disp(me.message);
                    r_param_s(hind,:) = r_param_s_copy;
                    r_param(hind,:) = r_param_copy;
                    r_error_s(hind,:) = r_error_s_copy;
                    r_error(hind,:) = r_error_copy;
                end

            end
            
            try
                % the molecular ions and neutral NO
                outputs = fc.get_point(glat,glon,galt,r_param(hind,1),r_param(hind,2)*r_param(hind,3),r_param(hind,2));
                O2p(hind) = outputs{5}/r_param(hind,1);
                NOp(hind) = outputs{6}/r_param(hind,1);
                NO(hind) = outputs{9};
                if isnan(O2p(hind)) | isnan(NOp(hind))
                    error('NaN composition, skipping the flipchem fit.');
                end
            catch me
                disp(me.message);
            end
            
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  smoohing in altitude ##
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        % at least 3 gates are needed for smoothing in range
        if nhei > 2
            
            % Form a correlation prior in range (height) direction
            A = zeros((nhei-1)+(nhei-2),nhei);        
            SNe = A(:,1);        
            STi = A(:,1);        
            STr = A(:,1);        
            SVi = A(:,1);        
            SOp = A(:,1);        
            SVi = A(:,1);
            
            Aind = 1;
            
            
            % The correlation powers solved from known variances, height steps, and correlation lengths
            corrP = r_error_s(:,1:6).^2;
            corrP(:,1) = corrP(:,1).*dheights./(hsNe*hsAlt);
            corrP(:,2) = corrP(:,2).*dheights./(hsTi*hsAlt);
            corrP(:,3) = corrP(:,3).*dheights./(hsTr*hsAlt);
            corrP(:,4) = NaN;
            corrP(:,5) = corrP(:,5).*dheights./(hsVi*hsAlt);
            corrP(:,6) = corrP(:,6).*dheights./(hsOp*hsAlt);
            
            % The first order terms.
            % M is always zero for the first and higher order terms
            % The zeroth-order terms are added later
            for hind = 1:(nhei-1)
                
                A(Aind,hind) = 1;
                
                A(Aind,hind+1) = -1;
                
                SNe(Aind) = 2 * corrP(hind,1) * dheights(hind) / (hsNe*hsAlt(hind));
                
                STi(Aind) = 2 * corrP(hind,2) * dheights(hind) / (hsTi*hsAlt(hind));
                
                STr(Aind) = 2 * corrP(hind,3) * dheights(hind) / (hsTr*hsAlt(hind));
                
                SVi(Aind) = 2 * corrP(hind,5) * dheights(hind) / (hsVi*hsAlt(hind));
                
                SOp(Aind) = 2 * corrP(hind,6) * dheights(hind) / (hsOp*hsAlt(hind));
                
                Aind = Aind+1;
            end
            
            % The second order terms
            % NOTE: This is approximately OK also when the altitude resolution changes, because we assume that
            % the parameters are constant within a gate...
            for hind = 2:(nhei-1)
                
                A(Aind,hind-1) = 1;
                A(Aind,hind) = -2;
                A(Aind,hind+1) = 1;
                
                SNe(Aind) = 8 * corrP(hind,1) * (dheights(hind) / (hsNe*hsAlt(hind)))^3;
                
                STi(Aind) = 8 * corrP(hind,2) * (dheights(hind) / (hsTi*hsAlt(hind)))^3;
                
                STr(Aind) = 8 * corrP(hind,3) * (dheights(hind) / (hsTr*hsAlt(hind)))^3;
                
                SVi(Aind) = 8 * corrP(hind,5) * (dheights(hind) / (hsVi*hsAlt(hind)))^3;
                
                SOp(Aind) = 8 * corrP(hind,6) * (dheights(hind) / (hsOp*hsAlt(hind)))^3;
                
                Aind = Aind+1;
            end
            
            
            % Combine all parametes in one large theory matrix
            [n1,n2] = size(A);
            
            % We have five parameters
            Acomb = zeros(5*n1,5*n2);
            Scomb = NaN(5*n1,1);
            for ipar = 1:5
                Acomb( ((ipar-1)*n1+1) : (ipar*n1) , ((ipar-1)*n2+1) : ...
                       (ipar*n2) ) = A;
            end
            Scomb = [ SNe ; STi ; STr ; SVi ; SOp ];
            Qcomb = Acomb'*diag(1./Scomb)*Acomb;

        else
            % zero information if we did not smooth in range
            Qcomb = zeros(nhei*5,nhei*5);
        end
            
        %
        % The zeroth order terms are the measurements from the previous step. These are with a
        % covariance matrix in each gate
        %
        Cfit = zeros(nhei*5,nhei*5);
        Mfit = NaN(nhei*5,1);
        for ihei = 1:nhei
            fitcov = vec2covm(r_error_s(ihei,:));
            Cfit( ((0:4)*nhei+ihei) , ((0:4)*nhei+ihei) ) = fitcov([1 2 3 5 6],[1 2 3 5 6]);
            Mfit( ((0:4)*nhei)+ihei ) = r_param_s(ihei,[1 2 3 5 6]);
        end
        
        % Precision matrix
        Qfit = inv(Cfit);

        % Solve the whole problem (zeroth + first + second order)
        Cpost = inv( Qfit + Qcomb);
        Xpost = Cpost * Qfit * Mfit;
        
        % Pick the parameter profiles
        NeCorr = Xpost(1:nhei);
        TiCorr = Xpost((nhei+1):(2*nhei));
        TrCorr = Xpost((2*nhei+1):(3*nhei));
        ViCorr = Xpost((3*nhei+1):(4*nhei));
        OpCorr = Xpost((4*nhei+1):(5*nhei));
        
        % Standard deviations. NOTICE: we could pick the full covariance matrices at each height!
        NeErrCorr = sqrt(diag(Cpost(1:nhei,1:nhei)));
        TiErrCorr = sqrt(diag(Cpost((nhei+1):(2*nhei),(nhei+1):(2*nhei))));
        TrErrCorr = sqrt(diag(Cpost((2*nhei+1):(3*nhei),(2*nhei+1):(3*nhei))));
        ViErrCorr = sqrt(diag(Cpost((3*nhei+1):(4*nhei),(3*nhei+1):(4*nhei))));
        OpErrCorr = sqrt(diag(Cpost((4*nhei+1):(5*nhei),(4*nhei+1):(5*nhei))));

        
        % Collect the smoothed parameters and convert back to physical units
        partmp = scaled_to_real([NeCorr(:),TiCorr(:),TrCorr(:),r_param_s(:,4),ViCorr(:),OpCorr(:)]);
        %errtmp = scaled_to_real([NeErrCorr(:),TiErrCorr(:),TrErrCorr(:),r_error_s(:,4),ViErrCorr(:),OpErrCorr(:)]);


        % pick the full covariance matrices and conver to physical units
        errtmp_s = r_error_s;
        % indices of fitted parameters in the full arrays
        ipar1 = [1 2 3 5 6];
        for ihei = 1:nhei
            % indices of fitted parameters from this gatesin Cpost
            ipar2 = (0:4)*nhei + ihei;
            % initialize the covariance matrix
            covarhei = vec2covm(r_error_s(ihei,:));
            % replace covariances of the fitted parameters with values from the smoothing
            covarhei(ipar1,ipar1) = Cpost(ipar2,ipar2);
            % convert back to GUISDAP vector format
            errtmp_s(ihei,:) = covm2vec(covarhei);
        end
        % back to physical units
        errtmp = scaled_to_real(errtmp_s);

        
        % skip the smoothing if the values are insane
        if ~any(any(isnan(partmp))) & ~any(any(isnan(errtmp))) & all(all(imag(partmp)==0)) & all(all(imag(errtmp)==0)) & ~any(any(errtmp(:,1:6)<0))
            
            % Copy the smoothed profiles to the final prior matrix
            apriori2(:,1) = partmp(:,1);
            apriori2(:,2) = partmp(:,2);
            apriori2(:,3) = partmp(:,3);
            apriori2(:,5) = partmp(:,5);
            apriori2(:,6) = partmp(:,6);
            
            % The final prior variances
            apriorierror2(:,1) = sqrt( errtmp(:,1).^2 + tsNe^2*tstep);
            apriorierror2(:,2) = sqrt( errtmp(:,2).^2 + tsTi^2*tstep);
            apriorierror2(:,3) = sqrt( errtmp(:,3).^2 + tsTr^2*tstep);
            apriorierror2(:,4) = 0;
            apriorierror2(:,5) = sqrt( errtmp(:,5).^2 + tsVi^2*tstep);
            apriorierror2(:,6) = sqrt( errtmp(:,6).^2 + tsOp^2*tstep);

        else
            % if something failed we skip the smoothing
            Cpost = Cfit;
            Qpost = Qfit;
            partmp = r_param;
            errtmp = r_error;

            % print a warning
            disp('Smoothing in altitude failed, skipping the smoothing.')
            
        end
        
        % boundary conditions
        for ihei = 1:nhei
            % Ne
            if heights(ihei)<hlimNe(1) | heights(ihei)>hlimNe(2)
                apriori2(ihei,1) = apriori(ihei,1);
                apriorierror2(ihei,1) = 1e-3;
            end
            % Ti
            if heights(ihei)<hlimTi(1) | heights(ihei)>hlimTi(2)
                apriori2(ihei,2) = apriori(ihei,2);
                apriorierror2(ihei,2) = 1e-3;
            end
            % Tr
            if heights(ihei)<hlimTr(1) | heights(ihei)>hlimTr(2)
                apriori2(ihei,3) = apriori(ihei,3);
                apriorierror2(ihei,3) = 1e-3;
            end
            % Vi
            if heights(ihei)<hlimVi(1) | heights(ihei)>hlimVi(2)
                apriori2(ihei,5) = apriori(ihei,5);
                apriorierror2(ihei,5) = 1e-3;
            end
            % O+
            if heights(ihei)<hlimOp(1) | heights(ihei)>hlimOp(2)
                apriori2(ihei,6) = apriori(ihei,6);
                apriorierror2(ihei,6) = 1e-3;
            end
        end

        % The matrix G in eq 8.6 of Särkkä. This is stored in the data files and used in the smoothing step
        scaledprederr = real_to_scaled(apriorierror2);
        predvars = scaledprederr(:,[1 2 3 5 6]).^2;
        BAFIM_G = Cfit * (Cpost * Qfit)' * diag( 1./predvars(:));

        % copy the parameters smoothed in altitude in matrices that are appended to the data files
        r_param_rcorr = r_param;
        r_param_rcorr(:,1:6) = partmp(:,1:6);
        %        r_error_rcorr = r_error;
        %        r_error_rcorr(:,1:6) = errtmp(:,1:6);
        r_error_rcorr = errtmp;
        r_param_filter = r_param;
        r_error_filter = r_error;

        r_dp = r_param(:,6);

        % write the matrix G and the smoothed parameters in the GUISDAP output file.
        % Write also r_param and r_error, since they have changed in the flipchem fit.
        % Write an identical copy of r_param and r_error in r_param_filter and r_error_filter.        
        outfile = fullfile(result_path,filename);

        savesuccess = false;
        itry = 0;
        while ~savesuccess
            savesuccess = true;
            try
                save(outfile,'BAFIM_G','r_param','r_error','r_status','r_param_filter','r_error_filter','r_param_rcorr','r_error_rcorr','r_dp','O2p','NOp','NO','kdeb','-append');
            catch
                savesuccess = false;
                itry = itry + 1;
                pause(5);
                if itry > 12
                    error(['cannot write to file ' outfile])
                end
            end
        end


        % merge the guisdap output files into one large file to avoid too many files...
        if merge_output_files
            merge_mat(result_path,true,true);
        end

    else
        % Prior for the very first time step or after a long data gap from the IRI model. apriori2 already contains the IRI parameters, except for Ne which might be from power profiles, just form the error array here.


        % Ne
        apriorierror2(:,1) = tsNe*sqrt(tstep)./hsAlt; % smaller process noise on the smooth and low Ne topside
        apriorierror2(heights<hlimNe(1),2) = 1e-3;
        apriorierror2(heights>hlimNe(2),2) = 1e-3;
        % Ti
        Apriorierror2(:,2) = tsTi*sqrt(tstep);
        apriorierror2(heights<hlimTi(1),2) = 1e-3;
         apriorierror2(heights>hlimTi(2),2) = 1e-3;
        % Tr
        apriorierror2(:,3) = tsTr*sqrt(tstep);
        apriorierror2(heights<hlimTr(1),3) = 1e-3;
        apriorierror2(heights>hlimTr(2),3) = 1e-3;
        % coll
        apriorierror2(:,4) = 0;
        % Vi
        apriorierror2(:,5) = tsVi*sqrt(tstep);
        apriorierror2(heights<hlimVi(1),5) = 1e-3;
        apriorierror2(heights>hlimVi(2),5) = 1e-3;
        % composition
        apriorierror2(:,6) = tsOp*sqrt(tstep);
        apriorierror2(heights<hlimOp(1),6) = 1e-3;
        apriorierror2(heights>hlimOp(2),6) = 1e-3;

        % if istep > 1 the previous output file should be deleted, because we want to completely skip that time step
        if istep > 1
            outfile = fullfile(result_path,filename);
            try
                disp(['Deleting' outfile])
                delete(outfile);
            catch
                disp('delete failed');
            end
        end
        
    end
    
    % make sure that the prior values are within our limits, this should rarely have any effect
    for ihei=1:nhei
        apriori2(ihei,1:6) = max(apriori2(ihei,1:6),paramlims(1,:));
        apriori2(ihei,1:6) = min(apriori2(ihei,1:6),paramlims(2,:));
    end

    
    % Copy the final prior model to aprioriprev and apriorierrorprev
    aprioriprev = apriori2;
    apriorierrorprev = apriorierror2;

    % update the output file name for the next integration period
    filename=sprintf('%08d.mat',fix(tosecs(d_time(2,:))));

    % increment the step counter
    istep = istep + 1;
    
    % Let the user know which prior is being used
    disp('apriorimodel_bafim_flipchem')
    
end

