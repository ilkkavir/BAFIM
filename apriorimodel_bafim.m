function [apriori2,apriorierror2] = apriorimodel_bafim(apriori,apriorierror,heights,fit_altitude)
%
%
% [apriori2,apriorierror2] = apriorimodel_bafim(apriori,apriorierror,heights)
%
% Bayesian filtering in time and correlation priors in range direction for Ne, Ti, Tr, Vi, and composition.
%
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
% See also: bafim_smoother
%
% IV 2018-2020
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
    
    global r_param r_error r_res r_status d_time path_GUP result_path v_Boltzmann v_amu
    
    % We need to pass some variables from one timestep to another
    persistent d_time_prev aprioriprev apriorierrorprev filename istep

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
    
    % The final correlation lengths are products hsXX*H, where H is the plasma scale height as calculated from IRI parameters.
    H = v_Boltzmann .* apriori(:,2).*(1+apriori(:,3))./2 ... 
        ./ ( v_amu .* ( apriori(:,6)*16 + 30.5*(1-apriori(:,6)) ) .* 9.82.*(6372./(6372+heights)).^2 );

    hsAlt = H/1000;
    
    % r_param is empty on the first iteration. Initialize d_time_prev, add the guisdap init directory to the path, and copy this file to the result directory.
    if size(r_param,1)==0
        d_time_prev = d_time(1,:);
        addpath([path_GUP 'init']);
        istep = 1;
    end
    if istep==2
        % copy the prior model in the output directory to help in development work
        % copy only on the second step to be sure that result_path has been properly parsed
        copyfile(which('apriorimodel_bafim.m'),fullfile(result_path,'apriorimodel_bafim.m'));
    end

    
    % Time step length
    tstep = seconds(datetime(d_time(2,:))-datetime(d_time_prev));
    d_time_prev = d_time(2,:);
    
    if nhei > 1
        % Approximate widths of the height gates
        dheights = diff(heights);
        dheights = [dheights(1) ; dheights];
    end
        
    % Copy the default prior model
    apriori2 = apriori;
    apriorierror2 = apriorierror;

    % if this is not the first time step
    if istep>1
        
        % Replace unreasonable values with the previous predictions. 
        for hind = 1:nhei
            if ~any(r_status(hind)==[0 3]) | any(isnan(r_param(hind,1:6))) | ...
                    r_res(hind,1)>100 | any(r_param(hind,1:5)<paramlims(1,1:5)) | ...
                    any(r_param(hind,1:5)>paramlims(2,1:5))

                r_param(hind,1:6) = aprioriprev(hind,1:6);
                r_error(hind,1:6) = apriorierrorprev(hind,1:6);
                r_error(hind,7:end) = 0;
            end
        end
        
        % Convert to scaled units
        r_param_s = real_to_scaled(r_param);
        r_error_s = real_to_scaled(r_error);
        apriori_s = real_to_scaled(apriori);
        apriorierror_s = real_to_scaled(apriorierror);

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
            Cfit( ((0:4)*nhei+ihei) , ((0:4)*nhei+ihei) ) = fitcov([1 2 ...
                                3 5 6],[1 2 3 5 6]);
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
        
        % Standard deviations
        NeErrCorr = sqrt(diag(Cpost(1:nhei,1:nhei)));
        TiErrCorr = sqrt(diag(Cpost((nhei+1):(2*nhei),(nhei+1):(2*nhei))));
        TrErrCorr = sqrt(diag(Cpost((2*nhei+1):(3*nhei),(2*nhei+1):(3*nhei))));
        ViErrCorr = sqrt(diag(Cpost((3*nhei+1):(4*nhei),(3*nhei+1):(4*nhei))));
        OpErrCorr = sqrt(diag(Cpost((4*nhei+1):(5*nhei),(4*nhei+1):(5*nhei))));
                
        % Convert the prediction back to physical units
        partmp = scaled_to_real([NeCorr(:),TiCorr(:),TrCorr(:),r_param_s(:,4),ViCorr(:),OpCorr(:)]);
        errtmp = scaled_to_real([NeErrCorr(:),TiErrCorr(:),TrErrCorr(:),r_error_s(:,4),ViErrCorr(:),OpErrCorr(:)]);
        
        % Copy the smoothed profiles to the final prior matrix
        apriori2(:,1) = partmp(:,1);
        apriori2(:,2) = partmp(:,2);
        apriori2(:,3) = partmp(:,3);
        apriori2(:,5) = partmp(:,5);
        apriori2(:,6) = partmp(:,6);
        
        % The final prior variances
        apriorierror2(:,1) = errtmp(:,1) + tsNe*sqrt(tstep);
        apriorierror2(:,2) = errtmp(:,2) + tsTi*sqrt(tstep);
        apriorierror2(:,3) = errtmp(:,3) + tsTr*sqrt(tstep);
        apriorierror2(:,4) = 0;
        apriorierror2(:,5) = errtmp(:,5) + tsVi*sqrt(tstep);
        apriorierror2(:,6) = errtmp(:,6) + tsOp*sqrt(tstep);

        
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

        % copy the smoothed parameters in matrices that are appended to the data files
        r_param_rcorr = r_param;
        r_param_rcorr(:,1:6) = partmp(:,1:6);
        r_error_rcorr = r_error;
        r_error_rcorr(:,1:6) = errtmp(:,1:6);
        r_param_filter = r_param;
        r_error_filter = r_error;

        % write the matrix G and the smoothed parameters in the GUISDAP output file.
        % Write also r_param and r_error, since they might have changed if there were failed fits.
        % Write an identical copy of r_param and r_error in r_param_filter and r_error_filter.
        outfile = fullfile(result_path,filename);
        save(outfile,'BAFIM_G','r_param','r_error','r_param_filter','r_error_filter','r_param_rcorr','r_error_rcorr','-append');
        
    else
        % Prior for the very first time step from the IRI model. apriori2 already contains the IRI parameters, just form the error array here. 
        
        % Ne
        apriorierror2(:,1) = tsNe*sqrt(tstep);
        apriorierror2(heights<hlimNe(1),2) = 1e-3;
        apriorierror2(heights>hlimNe(2),2) = 1e-3;
        % Ti
        apriorierror2(:,2) = tsTi*sqrt(tstep);
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
    disp('apriorimodel_bafim')
    
end

