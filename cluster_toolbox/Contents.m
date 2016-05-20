% Cluster toolbox, a collection of chemometrics functions 
% type help function_name for detailed instruction for each functions
% Copyright <a href="http://www.biospec.net">Laboratory for Bioanalytical Spectroscopy</a>
%
%Signal processing functions
%--------------------------------------------------------------------------
% asysm         asymmetrical least square baseline estimation (main
%               function) (new!)
% baseline_correction
%               interactive baseline correction usng asymmetrical least
%               squares.(new!)
% band_area     integrate a peak/spectral band interactively (new!)
% band_area2    integrate a peak/spectral band automatically by specifying 
%               its start and end points (new!)
% CO2corr       corrects for CO2 interference in FT-IR data
% derivats      Savitzky & Golay derivatives main function
% detrendm      detrends from first and last bin of a spectral matrix
% emsc          extended multiplicative scatter correction
% gaussian_smooth
%               Gaussian smooth function (new!)
% mos           morphological factors/scores, multivariate signal to noise
%               ratio measurements (new!)
% normal        min-max normalisation (to a range of 0 to 1000)
% normhigh      normalise to the highest point in each spectrum
% normtot       normalise to sum of total intensities
% dosc          direct orthogonal signal correction using Westerhuis et. al
%               algorithm (new!)
% osc           orthogonal signal correction using Fearn's algorithm (new!)
% qc_corr       correct long-term instrumental drifting using QCs (new!)
% scalem        min-max normalisation (to a range of 0 to 1)
% vecnorm       unit vector normalization (sum of squares of each row
%               vector/spectrum = 1)
%
%Modelling functions
%--------------------------------------------------------------------------
% asca          ANOVA-simultaneous component analysis with permutation
%               test for validation (new!)
% cpca          conensus principal component analysis, a multi-block PCA
%               model with emphasis on capturing major variance (new!)
% dfa           discriminant function analysis
% hpca          hierarchical principal component analysis, another
%               multi-block PCA with emphasis on capturing common trend
%               (new!)
% msca          multilevel simultaneous component analysis (new!)
% opls          orthognal projection to latent structure PLS 1 model (new!)
% opls2         orthognal projection to latent structure PLS 2 model (new!)
% oplsda        orthognal projection to latent structure for discriminant
%               anaysis (new!)
% pca           principal component analysis using NIPALS algorithm
% pls           partial least square analysis, support both PLS1 and PLS2
%               (new!)
% plspred       predict the outcome of new data using pls model, also
%               outputs the projected pls scores (new!)
% plspred2      same as plspred2 but faster, also the outputs regression
%               coefficients (new!)
% plsda_boots   partial least square for discriminant analysis with build
%               in bootstrapping validation (new!)
% plsr_boots    partial least sqaure regression with build in bootstrapping
%               validation (new!)
% projpcdf      pc-dfa using traing data then projection of test data into
%               the pc-dfa space
% oc_clustering performs hierarchical clustering using the OC program by 
%               Dr. Geoffrey J. Barton, this is a Matlab frontend of ocnt.exe
%
%Plot functions
%--------------------------------------------------------------------------
% error_ellipse plot an error ellipse, or ellipsoid, defining a given 
%               confidence region (new!)
% gradientclass_plot
%               2D plot with a gradient colour map, useful for data with a
%               ordinal changing conditions (e.g. time, temperature etc.) 
%               (new!) 
% gradientclass_plot2
%               2D plot for two types of classes (e.g. two different
%               experimental conditions), one differentiated by symbols and
%               another by different colours. (new!)
% multiclass_plot
%               use a random combination of symbols, colours to present a
%               large number of classes (up to 155 different classes)(new!)
% p2d_col       plots a 2D plot with a specific column(variable) against another
% plot_dfa      plot dfa scores with labels, titles, etc
% plot_map      plot a map (contour, surface or both) according to a matrix
% plot_pca      plot PCA with labels, title, etc
% plotftir      plot axis correctly for FT-IR in colour
% plotnm2       plots a 2D plot with numbers from 1 to number of rows
% plotnm3       plots a 3D plot with numbers from 1 to number of rows
% plotpyms      plots a single MS spectrum as % total
% pltnmred      same as plotnm2
%
% Data I/O
%--------------------------------------------------------------------------
% get_spc       Read Thermo-Galactic SPC format file into MATLAB (need
%               bioinformatics toolbox) (new!)
% get_cdf       Read low resolution GC/LC-MS netcdf format files into MATLAB 
%               (warning! m/z will be rounded to unit mass!) and output as 
%               a matrix (new!)
% get_cdf_tic   Read Total Ion Currents (TIC) from netcdf into MATLAB (new!)
% get_cdf_lcms  Read high resolution GC/LC-MS netcdf format files into
%               MATLAB in cell array without rounding.
% get_cdf_dims  Read high resolution direct-infusion MS netcdf files into
%               MATLAB in cell array without rounding.
% get_maldi_text
%               Read exported MALDI-T.o.F.-MS data in ASCII into MATLAB
%               (new!)
%
%Support functions (called by other functions, normally shouldn't be used alone)
%--------------------------------------------------------------------------
% dfa_red       same as plot_dfa, called by projpcdf
% difsm         support function for asysm, called by asysm
% difsmw        support function for asysm, called by asysm
% func_accept   support function for baseline_correction, called by
%               baseline_correction
% func_show     support function for baseline_correction, called by
%               baseline_correction    
% explv         calculate explained variance% of a pca model, called by PCA
% filtnam       sort the data according to a pre-defined experiment design
% floatout      outputs floating point data of 12.8
% gendst        general distance calculation procedure
% groups        returns a group files for dfa in order
% lintrans      make a linear transformation from one range to another
% makeidx       returns an index (idx) files base on find numbers you want in a groups file
% makemage      write a file in kinemage format
% meanidx       calculates means of a matrix based on index from a group file
% meanit        calculates means based on a pre-defined experiment design
% names         returns a names (integer character) files for plotting in order
% namesidx      orders names in order of groups on first occurence
% plot_df1      same as plotftir
% plot_df2      same as plotftir
% plotdetr      same as plotftir
% pca_red       same as plot_pca, called by projpcdf
% deriver       Savitzky & Golay derivatives engine, called by derivats
% tw_gen        DFA engine, called by dfa
% num2int       converts numerical matrix M into integer/character matrix
