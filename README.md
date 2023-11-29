Bayesian Filtering Module (BAFIM)


BAFIM is an additional module to the GUISDAP incoherent scatter analysis tool ( https://www.eiscat.se/scientist/user-documentation/guisdap-9-0/ ). BAFIM replaces the standard IRI-based prior model with Bayesian filtering in time. Correlation priors are used in the prediction step of the filter to keep the predicted parameter profiles smooth in range direction in field-aligned measurements. The module can be used also with beam directions oblique to B and in analysis of remote data.

Reference: Virtanen, I.I, Tesfaw, H, Roininen, L., Lasanen, S., and Aikio, A., Bayesian filtering in incoherent scatter plasma parameter fits, manuscript submitted to JGR space physics, 2020. 

Ilkka Virtanen (ilkka.i.virtanen@oulu.fi) 2020



#################################################################################################
INSTALLATION:

1. copy the files somewhere in your MATLAB path, for example the GUISDAP anal/ directory.

2. add the following lines in the end of the standard GUISDAP file anal/ionomodel.m:

   apriorimodel=['apriorimodel_' iono_model];
   if exist([apriorimodel '.m'])
       [apriori,apriorierror] = feval(apriorimodel,apriori,apriorierror,heights,fit_altitude);
   end



#################################################################################################
USAGE:

The scaling factors etc. mentioned below are explained in the reference paper. 

To use BAFIM and Flipchem in GUISDAP fits:

1. Select 'bafim_flipchem' as the ionospheric model in GUISDAP, i.e. write the following in the 'Special' box of the 'GUISDAP for dummies' window:

iono_model='bafim_flipchem'

2. Set the lowest (hmin) and highest (hmax) altitude where each parameter is fitted, and scaling factors for process noise (st) and correlation length (sh) using the first four colums and six rows of the 'fit_altitude' array:

fit_altitude(1:6,1:4) = [ hmin_Ne , hmax_Ne , sh_Ne , st_Ne;
		      	  hmin_Ti , hmax_Ti , sh_Ti , st_Ti;
			  hmin_Tr , hmax_Tr , sh_Tr , st_Tr;
			     0    ,    0    ,   0   ,   0  ;
			  hmin_Vi , hmax_Vi , sh_Vi , st_Vi;
			  hmin_Op , hmax_Op , sh_Op , st_Op ]

For example:

fit_altitude(1:6,1:4) = [   0 , Inf ,  0.05 , 2.5e11 ;
		      	   80 , Inf ,  0.1  ,     30 ;
			   97 , Inf ,  0.1  ,   0.03 ;
			    0 ,   0 ,    0  ,      0 ;
			   80 , Inf , 0.05  ,    2.5 ;
			  150 , 500 , 0.05  ,   0.01 ]

3. Tell the GUISDAP satellite check routine not to cut the profiles and set transmitter phase pushing to zero:

a_satch.cut=0
a_phasepush=0

3. Start GUISDAP analysis as usual (hit 'GO')




To run the RTS smoother, use the command:

   bafim_smoother(<guisdap-output-dir>)

After the smoothing step the Bayesian smoothing results are stored in variables r_param and r_error, as well as in r_param_smooth and r_error_smooth, while the Bayesian filtering outputs are in r_param_filter and r_error_filter. The results stored in the data files under r_param and r_error may be selected using the function bafim_select. 
