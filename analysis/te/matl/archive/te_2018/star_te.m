clear

str=computer;
if strcmp(str,'PCWIN64')
    %restoredefaultpath
    cd C:\Users\Dobri\Desktop
    addpath('C:\Users\Dobri\Documents\MATLAB\fieldtrip')
    ft_defaults
    addpath('C:\Users\Dobri\Documents\MATLAB\TRENTOOL3\')
    addpath('C:\Users\Dobri\Desktop\DGD\STAR_data\te')
    load('C:\Users\Dobri\Desktop\dgd_larger_data_files\head_bobbing_lv_data_for_spike_sync_surrogates_2017-09-06-16-40.mat','lvpos','lvtime','lvmeta','censure')
    OutputDataPath = pwd; %'results\';
else
    %restoredefaultpath
    cd ~/Desktop
    addpath('~/matlab-playwith/fieldtrip/')
    ft_defaults
    addpath('~/matlab-playwith/TRENTOOL3/')
    addpath('/home/dobri/mcmc/STAR_data/te');
    load('~/Desktop/head_bobbing_lv_data_for_spike_sync_surrogates_2017-09-06-16-40.mat','lvpos','lvtime','lvmeta','censure')
    OutputDataPath = pwd; %'results/';
end

for t=1:8
    lvpos{t}(:,:,censure)=[];
end
selected_lv=1;

%% define cfg for TEprepare.m
cfgTEP = [];

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 10;      % minimum u to be scanned
cfgTEP.predicttimemax_u    = 500;	  % maximum u to be scanned
cfgTEP.predicttimestepsize = 10; 	  % time steps between u's to be scanned

% estimator
cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% ACT estimation and constraints on allowed ACT(autocorelation time)
cfgTEP.actthrvalue = 1000;   % threshold for ACT
cfgTEP.maxlag      = 200; %?
% cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials
cfgTEP.trialselect = 'no';

% optimizing embedding
cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragdim         = 2:6;       % criterion dimension  %?
cfgTEP.ragtaurange    = [0.2 0.4]; % range for tau
cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps

% kernel-based TE estimation
cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
cfgTEP.sizeNei = 4;                 % neigbours to analyse


%% define cfg for TEsurrogatestats_ensemble.m
cfgTESS = [];

% use individual dimensions for embedding
cfgTESS.optdimusage = 'indivdim';

% statistical and shift testing
cfgTESS.tail           = 1;
cfgTESS.numpermutation = 5e4;
cfgTESS.shifttest      = 'no';
cfgTESS.shifttesttype  = 'TEshift>TE';
% cfgTESS.surrogatetype  = 'trialshuffling';
% cfgTESS.surrogatetype  = 'blockresampling'; %X % trialreverse
cfgTESS.surrogatetype  = 'blockreverse1'; % V
% cfgTESS.surrogatetype  = 'blockreverse2';
% cfgTESS.surrogatetype  = 'blockreverse3';

% results file name
cfgTESS.fileidout  = strcat(OutputDataPath,'stardata_');

TGA_results_trials=cell(1,8);
for trial=1:8
    %% convert data to the right format.
    data=[];
    data.trial{1}=squeeze(lvpos{trial}(:,selected_lv,:))';
    data.time{1}=lvtime{trial}';
    data.fsample=lvmeta.fr;
    for pp=1:size(lvpos{trial},3);data.label{1,pp}=['P' num2str(pp+(pp>=censure),'%02.f')];end
    %clear pp
    toi            = [min(data.time{1,1}),max(data.time{1,1})]; % time of interest
    repPred        = round(size(data.trial{1},2)*(1/2));    % size(data.trial{1,1},2)*(3/4);   %?
    channel        = data.label;  % channels to be analyzed
    
    %% Here we go!
    TGA_results = basic_TEpermtest(cfgTEP,cfgTESS,data,toi,repPred,channel);
    TGA_results_trials{trial} = TGA_results;
    % TGA_results = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);
    
    %% Save
    save(fullfile(OutputDataPath,['TGA_results' '_Trial' num2str(trial,'%02.f') '.mat']),'TGA_results');
    TEM=nan(33);k=0;for r=1:33;for c=1:33;if r~=c;k=k+1;TEM(r,c)=TGA_results.TEmat(k);end;end;end;imagesc(TEM,[-.15 -.01]);axis square
    print('-djpeg','-r300',fullfile(OutputDataPath,['TEM' '_Trial' num2str(trial,'%02.f') '.jpeg']))
    sTEM=nan(33);k=0;for r=1:33;for c=1:33;if r~=c;k=k+1;sTEM(r,c)=TGA_results.TEpermvalues(k,2);end;end;end;imagesc(sTEM);axis square
    print('-djpeg','-r300',fullfile(OutputDataPath,['STEM' '_Trial' num2str(trial,'%02.f') '.jpeg']))
end

return



%% optional: perform a post hoc correction for cascade effects and simple common drive effects
cfgGA = [];
cfgGA.threshold = 3;
cfgGA.cmc       = 1;
TGA_results_GA = TEgraphanalysis(cfgGA,TGA_results);
save([OutputDataPath 'Lorenz_1-2_TGA_results_analyzed_GA.mat'],'TGA_results_GA');



%% plotting
cfgPLOT = [];
cfgPLOT.layout        = lay_Lorenz; 		% see fieldtrip's ft_prepare_layout.m
cfgPLOT.electrodes    = 'highlights';
cfgPLOT.statstype     = 1;   		% 1: corrected; 2:uncorrected; 3: 1-pval; 4:rawdistance
cfgPLOT.alpha         = 0.05;
cfgPLOT.arrowpos      = 1;
cfgPLOT.showlabels    = 'yes';
cfgPLOT.electrodes    = 'on';
cfgPLOT.hlmarker      = 'o';
cfgPLOT.hlcolor       = [0 0 0];
cfgPLOT.hlmarkersize  = 4;
cfgPLOT.arrowcolorpos = [1 0 0];

figure;
TEplot2D(cfgPLOT,TGA_results_GA)
