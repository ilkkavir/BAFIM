# Bayesian Filtering Module (BAFIM)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4033903.svg)](https://doi.org/10.5281/zenodo.4033903)

BAFIM is an additional module to the [GUISDAP](https://www.eiscat.se/scientist/user-documentation/guisdap-9-0/) incoherent scatter analysis tool. BAFIM replaces the standard IRI-based prior model with Bayesian filtering in time. Correlation priors are used in the prediction step of the filter to keep the predicted parameter profiles smooth in range direction in field-aligned measurements. The module can be used also with beam directions oblique to B and in analysis of remote data.

The [Flipchem](https://github.com/amisr/flipchem) [![DOI](https://zenodo.org/badge/DOI/10.5281/10.5281/zenodo.3688853.svg)](https://doi.org/10.5281/zenodo.3688853) model can be used to support F region O+ ion fraction fits. 




Ilkka Virtanen (ilkka.i.virtanen@oulu.fi) 2024


## INSTALLATION

1. Copy the files somewhere in your MATLAB path, for example the GUISDAP anal/ directory.

2. Add the following lines in the end of the standard GUISDAP function ionomodel in the file anal/ionomodel.m (before the last end statement)

```
apriorimodel=['apriorimodel_' iono_model];
if exist([apriorimodel '.m'])
  [apriori,apriorierror] = feval(apriorimodel,apriori,apriorierror,heights,fit_altitude);
end
```


## How to use the module?

There are two version of BAFIM available, the original BAFIM and the new BAFIM_flipchem. These are included as separate functions. 


**Notice that user inputs are scaled differently in the present versions of BAFIM and BAFIM_flipchem!**


### BAFIM (the original versoin)

The scaling factors etc. mentioned below are explained in the reference paper [Virtanen et al., 2021](https://doi.org/10.1029/2020JA028700). 

To use BAFIM in GUISDAP fits:

1. Select 'bafim' as the ionospheric model in GUISDAP, i.e. write the following in the 'Special' box of the 'GUISDAP for dummies' window:

    iono_model='bafim'

2. Set the lowest (hmin) and highest (hmax) altitude where each parameter is fitted, and scaling factors for process noise (st) and correlation length (sh) using the first four colums and six rows of the 'fit_altitude' array:

```
fit_altitude(1:6,1:4) = [ hmin_Ne , hmax_Ne , sh_Ne , st_Ne;  
                          hmin_Ti , hmax_Ti , sh_Ti , st_Ti;  
                          hmin_Tr , hmax_Tr , sh_Tr , st_Tr;  
                             0    ,    0    ,   0   ,   0  ;  
                          hmin_Vi , hmax_Vi , sh_Vi , st_Vi;  
                          hmin_Op , hmax_Op , sh_Op , st_Op ]  
```

For example:
```
fit_altitude(1:6,1:4) = [   0 , Inf , 0.1 ,  2e11 ;  
                           80 , Inf , 0.3 ,    10 ;  
                          103 , Inf , 0.3 ,  0.05 ;  
                            0 ,   0 ,   0 ,     0 ;  
                           80 , Inf , 0.2 ,   2.5 ;  
                          130 , 350 , 0.2 , 0.003 ]  
```

3. Tell the GUISDAP satellite check routine not to create gaps in the lag profiles and set transmitter phase pushing to zero:

```
a_phasepush=0  
a_satch.cut=0
```

4. Start GUISDAP analysis as usual (hit 'GO')




To run the RTS smoother, use the command:

```
bafim_smoother(<guisdap-output-dir>)
```

After the smoothing step the Bayesian smoothing results are stored in variables r_param and r_error, as well as in r_param_smooth and r_error_smooth, while the Bayesian filtering outputs are in r_param_filter and r_error_filter. The results stored in the data files under r_param and r_error may be selected using the function bafim_select. 




### BAFIM_flipchem

**The Flipchem package must be installed**. BAFIM_flipchem calls the Flipchem python module. 

The scaling factors etc. mentioned below are explained in the reference paper [Virtanen et al., 2024](https://doi.org/10.22541/essoar.169945113.37163064/v1). **Notice that these are different from the original version** 

To use BAFIM_flipchem in GUISDAP fits:

1. Select 'bafim' as the ionospheric model in GUISDAP, i.e. write the following in the 'Special' box of the 'GUISDAP for dummies' window:

```
iono_model='bafim_flipchem'
```

2. Set the lowest (hmin) and highest (hmax) altitude where each parameter is fitted, and scaling factors for process noise (st) and correlation length (sh) using the first four colums and six rows of the 'fit_altitude' array:

```
fit_altitude(1:6,1:4) = [ hmin_Ne , hmax_Ne , sh_Ne , st_Ne;  
                          hmin_Ti , hmax_Ti , sh_Ti , st_Ti;  
                          hmin_Tr , hmax_Tr , sh_Tr , st_Tr;  
                             0    ,    0    ,   0   ,   0  ;  
                          hmin_Vi , hmax_Vi , sh_Vi , st_Vi;  
                          hmin_Op , hmax_Op , sh_Op , st_Op ]  
```

For example:

```
fit_altitude(1:6,1:4) = [   0 , Inf , 0.05 ,  2.5e11 ;  
                           80 , Inf , 0.10 ,      30 ;  
                   	   97 , Inf , 0.10 ,    0.03 ;  
                            0 ,   0 ,    0 ,       0 ;  
                      	   80 , Inf , 0.05 ,     2.5 ;  
                          150 , 500 , 0.05 ,    0.01 ]  
```

**Notice that the numbers are different from those in the orignal BAFIM, because the length scales scale with time step duration and the chemistry model supports the composition fits.**

3. Tell the GUISDAP satellite check routine not to create gaps in the lag profiles and set transmitter phase pushing to zero:

```
a_phasepush=0
a_satch.cut=0
```

4. Start GUISDAP analysis as usual (hit 'GO')



To run the RTS smoother, use the command:

```
bafim_flipchem_smoother(<guisdap-output-dir>)
```

Four different results are stored in each data file after the smoothing step,

1. r_param_filter and r_error_filter contain the Bayesian filtering outputs, i.e. the parameters immediately after GUISDAP fit and Flipchem correction. 
2. r_param_rcorr and r_error_rcorr contain the plasma parameters after applying the smoothing in range by correlation priors. 
3. r_param_smooth and r_error_smooth are Bayesian smoothing results calculated from r_param_filter and r_error_filter. These results were shown in [Virtanen et al., 2021](https://doi.org/10.1029/2020JA028700). 
4. r_param_rcorr_smooth and r_error_rcorr_smooth are Baeysian smoothing results calculated from r_param_rcorr and r_error_rcorr. These results are shown in [Virtanen et al., 2024](https://doi.org/10.22541/essoar.169945113.37163064/v1)

r_param_smooth and r_error smooth are written to r_param and r_error. However, this can be changed with the function bafim_flipchem_select to make any of these data easily readable by existing routines. 

## References

Ilkka I. Virtanen, Habtamu W. Tesfaw, Anita T. Aikio, Roger Varney, Antti Kero, and Neethal Thomas, F1 region ion composition in Svalbard during the International Polar Year 2007-2008, Submitted to Journal of Geophysical Research: Space Physics, 2023, [https://doi.org/10.22541/essoar.169945113.37163064/v1](https://doi.org/10.22541/essoar.169945113.37163064/v1)

Ilkka I. Virtanen, Habtamu W. Tesfaw, Lassi Roininen, Sari Lasanen, and Anita Aikio, Bayesian filtering in incoherent scatter plasma parameter fits, Journal of Geophysical Research: Space Physics, 126 (3), 2021, [https://doi.org/10.1029/2020JA028700](https://doi.org/10.1029/2020JA028700)

Richards, P.G., Bilitza, D., Voglozin, D. Ion density calculator (IDC): A new eﬃcient model of ionospheric ion densities. Radio Science, 45:1–11, 2010, [https://doi.org/10.1029/2009RS004332](https://doi.org/10.1029/2009RS004332).

Richards, P. G., Reexamination of ionospheric photochemistry. Journal of Geophysical Research: Space Physics, 116 (8), 1–15, [https://doi.org/10.1029/2011JA016613](https://doi.org/10.1029/2011JA016613)