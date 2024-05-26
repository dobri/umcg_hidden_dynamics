clear

%% set software paths
addpath('~/logos/umcg_hidden_dynamics/analysis/te/matl/')
addpath('~/logos/umcg_hidden_dynamics/analysis/te/matl/libs/TRENTOOL3')
addpath('~/logos/umcg_hidden_dynamics/analysis/te/matl/libs/fieldtrip');
ft_defaults;

%% define data paths
base_folder = '~/biomech/projects/side_projects/umcg/data/';
% base_folder = 'C:\Users\ddotov\biomech\projects\side_projects\umcg\data\';
OutputDataPath = fullfile(base_folder,'results/');
InputDataPath  = fullfile(base_folder,'raw_data/');

cd(fullfile(base_folder))
data_aggregated_already = 1;
if data_aggregated_already
    load(fullfile(InputDataPath,'DATA_2024-05-24.mat'))
else
    DATA = load_raw_in_DATA_per_task(InputDataPath);
end


%%
TGA_results = cell(1,numel(DATA));
%parfor task = 1:numel(TGA_results)
for task = 1:numel(TGA_results)
    data = DATA{task};

    %% define cfg for TEprepare.m
    cfgTEP = [];

    cfgTEP.toi                 = [min(DATA{task}.time{1,1}),max(DATA{task}.time{1,1})]; % time of interest
    cfgTEP.channel             = DATA{task}.label;

    % scanning of interaction delays u
    cfgTEP.predicttimemin_u    = 0;      % minimum u to be scanned
    cfgTEP.predicttimemax_u    = 1000;	  % maximum u to be scanned
    cfgTEP.predicttimestepsize = 200; 	  % time steps between u's to be scanned

    % estimator
    cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

    % ACT estimation and constraints on allowed ACT(autocorelation time)
    cfgTEP.trialselect = 'no';
    cfgTEP.actthrvalue = 10;   % threshold for ACT
    cfgTEP.maxlag      = 10;
    cfgTEP.minnrtrials = 1;   % minimum acceptable number of trials

    % optimizing embedding
    cfgTEP.optimizemethod = 'ragwitz';  % criterion used
    cfgTEP.ragdim         = DATA{task}.dim(1,1):DATA{task}.dim(1,1)*3;       % criterion dimension
    cfgTEP.ragtaurange    = [.25 .5];
    cfgTEP.ragtausteps    = 5;         % steps for ragwitz tau steps
    cfgTEP.repPred        = 1e3; % size(data.trial{1,1},2)*(1/2);
    % cfgTESS.dim = DATA{task}.dim(1,1);
    % cfgTESS.tau = DATA{task}.tau(1,1);

    % kernel-based TE estimation
    cfgTEP.flagNei = 'Mass';            % neigbour analyse type
    cfgTEP.sizeNei = 4;                 % neigbours to analyse

    % set the level of verbosity of console outputs
    cfgTEP.verbosity = 'info_minor';
    %cfgTEP.verbosity = 'none';

    % define cfg for TE
    cfgTESS = [];

    % use individual dimensions for embedding
    cfgTESS.optdimusage = 'maxdim';
    % cfgTESS.optdimusage = 'indivdim';

    % statistical and shift testing
    cfgTESS.tail           = 1;
    % cfgTESS.numpermutation = 'findDelay';
    cfgTESS.numpermutation = 1e3;
    cfgTESS.shifttest      = 'no';

    cfgTESS.surrogatetype = 'trialshuffling'; % 'trialreverse','blockreverse2','swapneighbors';

    % don't calculate MI additionally to TE
    cfgTESS.MIcalc = 0;

    % results file name
    cfgTESS.fileidout = strcat(OutputDataPath,'out_');

    %% TE compuation, scanning over specified values of u
    TGA_results{task} = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);


    %% Null distribution from surrogate pairs
    % # repetitions of re-analyzing the full data w/ shifted (surrogate) pairs.
    TGA_results{task}.ntrials = size(TGA_results{task}.TEmat,2);
    TGA_results{task}.trial_shuffle_surr_n_per_trial = 3;
    TGA_results{task}.TEmat_sur2 = nan(numel(cfgTEP.channel),0);

    for repeats = 1:TGA_results{task}.trial_shuffle_surr_n_per_trial
        trials2 = randperm(TGA_results{task}.ntrials);
        data2 = data;
        for tr = 1:numel(trials2)
            data2.trial{tr}(2,:) = data.trial{trials2(tr)}(2,:);
        end

        cfgTEP2 = cfgTEP;
        cfgTEP2.predicttimemin_u    = min(TGA_results{task}.cfg.predicttime_u);      % minimum u to be scanned
        cfgTEP2.predicttimemax_u    = max(TGA_results{task}.cfg.predicttime_u);	  % maximum u to be scanned
        if cfgTEP2.predicttimemin_u == cfgTEP2.predicttimemax_u
            cfgTEP2.predicttimemax_u = cfgTEP2.predicttimemin_u + 100;
        end
        cfgTEP2.predicttimestepsize = (cfgTEP2.predicttimemax_u - cfgTEP2.predicttimemin_u); 	  % time steps between u's to be scanned
        cfgTEP2.ragdim = max(TGA_results{task}.cfg.dim);

        TGAtemp = InteractionDelayReconstruction_calculate(cfgTEP2,cfgTESS,data2);
        TGA_results{task}.TGA_surrogate_n{repeats} = TGAtemp;
        TGA_results{task}.TEmat_sur2 = [TGA_results{task}.TEmat_sur2 TGAtemp.TEmat];
        fprintf('%6.2f%% surrogates done.\n',repeats/TGA_results{task}.trial_shuffle_surr_n_per_trial*1e2)
    end
end

for task = 1:numel(TGA_results)
    m = mean(TGA_results{task}.TEmat_sur2,2);
    sd = std(TGA_results{task}.TEmat_sur2,[],2);
    TGA_results{task}.TEmat_rescaled = (TGA_results{task}.TEmat - m)./sd;

    colorvec = hsv(2);
    for d = 1:2
        subplot(4,2,task)
        [n,edges] = histcounts(TGA_results{task}.TEmat(d,:)','Normalization','probability');
        [n0,edges0] = histcounts(TGA_results{task}.TEmat_sur2(d,:)','Normalization','probability');
        plot(edges(2:end)+diff(edges(1:2))/2,n,'-','Color',colorvec(d,:).*.8,'linewidth',1.5)
        hold on
        plot(edges0(2:end)+diff(edges0(1:2))/2,n0,'--','linewidth',2,'Color',colorvec(d,:))
    end
    title(['Task ' num2str(DATA{task}.task{1})])
    hold off
    legend('TE_{Stim->Human}','TE_{Surrogate,Stim->Human}','TE_{Human->Stim}','TE_{Surrogate,Human->Stim}')
    legend('boxoff')

    DATA{task}.te = TGA_results{task}.TEmat_rescaled;
end
% f = fullfile(OutputDataPath,['PDF_TE_rescaled_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
% print('-djpeg','-r300',f)


T = table('Size',[3000,14],'VariableTypes',{'string','int8','int8','string',...
    'int8','string','string','int8',...
    'double','double','double','double','double','double'});
T.Properties.VariableNames = {'FileName','PP','Day','Date',...
    'TaskCode','TrainingPhase','TrainingPhaseLabel','Visual',...
    'Score','C','tau','PitchError','TE12rescaled','TE21rescaled'};
counter = 0;
for task = 1:numel(TGA_results)
    for tr = 1:TGA_results{task}.ntrials
        counter = counter + 1;
        T{counter,1} = DATA{task}.trial_name{1,tr};
        T{counter,2} = DATA{task}.pp{1,tr};
        T{counter,3} = DATA{task}.day{1,tr};
        T{counter,4} = DATA{task}.date{1,tr};
        T{counter,5} = DATA{task}.task{1,tr};
        T{counter,6} = DATA{task}.training_phase{1,tr};
        T{counter,7} = DATA{task}.training_phase_label{1,tr};
        T{counter,8} = DATA{task}.visual{1,tr};
        T{counter,9} = DATA{task}.scores(1,tr);
        T{counter,10} = DATA{task}.scores(2,tr);
        T{counter,11} = DATA{task}.scores(3,tr);
        T{counter,12} = DATA{task}.scores(4,tr);
        T{counter,13} = TGA_results{task}.TEmat_rescaled(1,tr);
        T{counter,14} = TGA_results{task}.TEmat_rescaled(2,tr);
    end
end
T(counter+1:end,:) = [];


%{
save(fullfile(OutputDataPath,['out_TGA_results_' char(datetime('now','Format','yyyy-MM-dd')) '.mat']),'TGA_results','T','OutputDataPath','InputDataPath');
writetable(T,fullfile(OutputDataPath,['Scores_TE_' char(datetime('now','Format','yyyy-MM-dd')) '.csv']))
%}
