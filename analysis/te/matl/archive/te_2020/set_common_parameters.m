function [cfgTEP,cfgTESS] = set_common_parameters()

%% define cfg for TEprepare.m
cfgTEP = [];
cfgTEP.verbosity = 'none';%'info_major';

% data
cfgTEP.channel             = {'Stim','User'};  % channels to be analyzed

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 20;      % minimum u to be scanned
cfgTEP.predicttimemax_u    = 520;	  % maximum u to be scanned
cfgTEP.predicttimestepsize = 50; 	  % time steps between u's to be scanned

% estimator
cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% ACT estimation and constraints on allowed ACT(autocorelation time)     <-
% The faughtdoodles keep switching unit convention between ms and sample.
cfgTEP.actthrvalue = 100;   % threshold for ACT
cfgTEP.maxlag      = 100;
cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials
cfgTEP.trialselect = 'no';

cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragtaurange    = [0.2 0.5]; % range for tau
cfgTEP.ragtausteps    = 5;         % steps for ragwitz tau steps

% kernel-based TE estimation
cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
cfgTEP.sizeNei = 4;                 % neigbours to analyse

cfgTEP.ensemblemethod = 'no';

cfgTEP.TEparallel.parON   = 'no'; %'yes'; % for paralell computing tool.
cfgTEP.TEparallel.workers =   16; % number of workers for parallel computing.


%% define cfg
cfgTESS = [];

% use individual dimensions for embedding
cfgTESS.optdimusage = 'indivdim';

% statistical and shift testing
cfgTESS.tail           = 1;
cfgTESS.numpermutation = 0;

cfgTESS.shifttest      = 'no';
cfgTESS.shifttesttype  = 'TEshift>TE';

% Other surrogate types. These work with single surrogate and permuation test.
% cfgTESS.surrogatetype  = 'trialshuffling';
% cfgTESS.surrogatetype  = 'blockresampling';
% cfgTESS.surrogatetype  = 'blockreverse1';
% cfgTESS.surrogatetype  = 'blockreverse2';
% cfgTESS.surrogatetype  = 'blockreverse3';