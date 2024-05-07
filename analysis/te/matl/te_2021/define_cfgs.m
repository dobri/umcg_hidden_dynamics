function [cfgTEP,cfgTESS] = define_cfgs(data)

%% define cfg for TEprepare.m
cfgTEP = [];

cfgTEP.toi                 = [min(data.time{1,1}),max(data.time{1,1})]; % time of interest
%cfgTEP.sgncmb              = {'A1' 'A2';'A2' 'A1'};  % channels to be analyzed
cfgTEP.channel             = data.label;

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 10;      % minimum u to be scanned
cfgTEP.predicttimemax_u    = 210;	  % maximum u to be scanned
cfgTEP.predicttimestepsize = 20; 	  % time steps between u's to be scanned

% estimator
cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% ACT estimation and constraints on allowed ACT(autocorelation time)
cfgTEP.actthrvalue = 100;   % threshold for ACT
cfgTEP.maxlag      = 100;
cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials

% optimizing embedding
cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragdim         = 3:9;       % criterion dimension
%cfgTEP.ragdim         = 3;       % data.dim(1,1) criterion dimension
cfgTEP.ragtaurange    = [0.75 1.]; % range for tau
cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps
cfgTEP.repPred        = 100; % round(size(data.trial{1,1},2)*(1/2));      % size(data.trial{1,1},2)*(3/4);

% kernel-based TE estimation
cfgTEP.flagNei = 'Mass';           % neigbour analyse type
cfgTEP.sizeNei = 4;                 % neigbours to analyse

% set the level of verbosity of console outputs
cfgTEP.verbosity = 'info_minor';
cfgTEP.verbosity = 'none';

%% define cfg for TE
cfgTESS = [];

% use individual dimensions for embedding
% cfgTESS.optdimusage = 'maxdim';
cfgTESS.optdimusage = 'indivdim';

%cfgTESS.dim = data.dim(1,1);
%cfgTESS.tau = data.tau(1,1);
% cfgTESS.dim and cfgTESS.tau

% statistical and shift testing
cfgTESS.tail           = 1;
%cfgTESS.numpermutation = 1e3;
cfgTESS.numpermutation = 5e3;
cfgTESS.shifttest      = 'no';
cfgTESS.alpha          = 1e-3;
%cfgTESS.shifttesttype  ='TEshift>TE';

%cfgTESS.surrogatetype  = 'trialshuffling';
cfgTESS.surrogatetype  = 'trialreverse';
%cfgTESS.surrogatetype  = 'blockreverse2';
%cfgTESS.surrogatetype  = 'swapneighbors';

% don't calculate MI additionally to TE
cfgTESS.MIcalc = 0;

