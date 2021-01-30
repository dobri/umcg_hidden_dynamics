function DATA = convert_csv_to_ft_and_te(InputDataPath)
%%
%{
Loop through participant folders and trials. The process in batches.
Group the trials in a structure with levels for participant and stimulus 
type. Thus, batches of trials are being passed to the analysis.
Some analysis parameters will be constrained to be the same for all
trials in the current batch, that is with the same stimulus type and participant.
Also useful for parallel processing which is needed because the analysis takes a lot of time.
By stimulus type we mean which dynamical system generates the stimulus.
Instead of wasting computation time and allowing for error, we will set
the embedding dimension based on our prior knowledge of the task dynamic.
m=2 for pure oscillaton (phase oscillator) and m=3 for nonlinear, chaotic stim.
The raw data is also pre-analyzed to get tau (for the embedding delay).
Then, after doing all the heavy (expect hours!) computation, the raw 
results are saved in one large structure for book-keeping. The outcome
variable needs to be re-arranged in an array allowing further stats.
%}


%% Prepare the data, arrange in the right format, get the conditions for the file name, etc..

% Assuming that each particpant's data is in a separate folder that is 
% labeled in a consistent manner. Here for demonstration I just used data
% where I was the pilot and a copy of it (for testing parfor and
% consistency of surrogate testing).
folderCell = dir(fullfile(InputDataPath,'pp*'));

DATA = cell(1); % The number of cells in DATA will grow.
% This approach might lead to excessive memory usage if you try to stuff
% all participants and all trials in the same structure. How many
% participants have you got? It is not necessary to have all data in the
% same structure. You can, for example, redesign this loop and the one
% below for the core TE analysis by loading and analyzing participant per
% participant. Still, because this analysis takes FOREVER, you want to make
% use of as many parallel workers as you can get from your computer. So if 
% you have 4 logical CPUs you want to at least have DATA with 4 full cells 
% to analyze at a time.

counter = 0;
for f=1:length(folderCell)
    S = readtable(fullfile(folderCell(f).folder,folderCell(f).name,'scores'));
    f_trials = dir(fullfile(folderCell(f).folder,folderCell(f).name,'trial*'));
    
    % How many different stimulus types in the given folder?
    tasks_vec = [];
    for tr = 1:numel(f_trials)
        tasks_vec(tr) = str2double(f_trials(tr).name(29:30));
    end
    tasks = unique(tasks_vec);
    
    % Gather by stimulus type and then trials that are of the given type.
    for ta = 1:numel(tasks) % That's not needed here. Tasks are fully b/w.
        counter = counter + 1;
        counter_d = 0;
        for tr = 1:numel(f_trials)
            if tasks_vec(tr) == tasks(ta)
                trial_file_name = f_trials(tr).name;
                
                %% Trial conditions and parameters
                day = str2double(f_trials(tr).name(11:16));
                hour = str2double(f_trials(tr).name(18:23));
                cond=str2double(trial_file_name(regexpi(trial_file_name,'task')+4:regexpi(trial_file_name,'task')+6));
                auditory=str2double(trial_file_name(regexpi(trial_file_name,'aud')+3));
                visual=str2double(trial_file_name(regexpi(trial_file_name,'vis')+3));
                epsilon=str2double(trial_file_name(regexpi(trial_file_name,'eps')+3:regexpi(trial_file_name,'eps')+5))/100;
                cond_vec = [day hour cond auditory visual epsilon];


                %% Don't search but decide the embedding dimension: 2 or 3,
                % depending on the task space (oscillator or nonlinear sys)
                % Important! How are the conditions coded in the file names?
                % When I last looked, the code was such that everything > 20 is a complex stimulus.
                % 10s for a phase oscillator, 20s for Chua, 30s for Lorenz.
                % In the earlier pilot trials, however, what I'm using here
                % as an example, 3 was the argument for Chua. Check the 
                % arguments of the python program used in the experiment.
                if cond >=20 || cond==3
                    m_stim=3;
                    m_pp=3;
                else
                    m_stim=2;
                    m_pp=2;
                end
                
                
                %% Load the data from one trial.
                x = dlmread(fullfile(f_trials(tr).folder,f_trials(tr).name),',',1,0);
                x = x(x(:,1)>10,:); % Crop the first 10 seconds of the trial!
                t = x(:,1); % time
                t = t - min(t); % for simplicity, shift back to zero.
                if auditory == 1
                    x_stim = x(:,11); % the stimulus
                    x_user = x(:,12); % the participant
                else % assume it's a visual modality condition.
                    x_stim = x(:,16); % the stimulus
                    x_user = x(:,17); % the participant
                end
                
                fprintf('%8.2f,',cond_vec,max(t),1/mean(diff(t)));fprintf('\n')
                if max(t)<40;keyboard;
                %else;pause
                end

                
                %% A little bit of pre-processing: smooth and detrend.
                x_user = smooth(x_user,10,'moving');
                % linear detrend the participant movement cause this can 
                % mess up several of the parameters (inf corr time).
                b = [ones(size(x_user)) t]\x_user;
                x_user = x_user - b(2)*t;

                stimstd = std(x_stim);
                stimmean = nanmean(x_stim);
                x_stim = (x_stim-stimmean)./stimstd;
                x_user = (x_user-stimmean)./stimstd;

                
                %% Estimate tau even though for now we do not use it. :/
                % On theory the optimal tau in PSR should be about a
                % quarter of the average oscillation period for signals
                % that are mostly oscillatory (narrow band). This tau will 
                % make the original signal and its delayed copy orthogonal.
                % So, you don't even need to search for the best tau using
                % computationally intensive methods, you just pass this value.
                % I have the quarter period tau here but I don't actually
                % use it. It's a good idea to find a way to pass it to 
                % TRENTOOL, skipping default methods and savings time.
                tau_stim = tau_quarter_period(x_stim,t);
                tau_pp = tau_quarter_period(x_user,t);
                
                counter_d = counter_d + 1;
                DATA{counter}.scores(counter_d,:) = S{tr,2:end};
                DATA{counter}.pp{counter_d} = folderCell(f).name;
                DATA{counter}.file_name{counter_d} = trial_file_name;
                DATA{counter}.conditions(counter_d,:) = cond_vec;
                DATA{counter}.trial{counter_d} = [x_stim x_user]';
                DATA{counter}.time{counter_d} = t';
                DATA{counter}.fsample(counter_d) = 1/mean(diff(t));
                DATA{counter}.label = {'Stim','User'};
                DATA{counter}.tau(:,counter_d) = [tau_stim;tau_pp];
                DATA{counter}.dim(:,counter_d) = [m_stim;m_pp];
            end
        end
        
        %% Just in case crop trials to same length. That's how fieldtrip wants them.
        % !!! Careful if you happen to have a short, interrupted trial in
        % the list! This will screw up everything.
        maxl = Inf;
        for tr = 1:size(DATA{counter}.time,2)
            maxl = min([maxl size(DATA{counter}.time{tr},2)]);
        end
        for tr = 1:size(DATA{counter}.time,2)
            DATA{counter}.trial{tr} = DATA{counter}.trial{tr}(:, 1:maxl);
            DATA{counter}.time{tr} = DATA{counter}.time{tr}(:, 1:maxl);
        end
    end
    fprintf('\n\n\n')
    %pause
    DATA{counter}.fsample = round(mean(DATA{counter}.fsample(counter_d)));
end

return

%% The core TE loop is going over batches of trials. Use the parfor on a multi-core CPU.
TE_rez = cell(1);
num_batches = numel(DATA);
parfor f = 1:num_batches
%for f=1%:num_batches
    fprintf('Analyzing files in %s, task %2.0f.\n',DATA{f}.pp{1},DATA{f}.conditions(1,3))

    % Load
    data=DATA{f};
    data.fsample=mean(data.fsample);
    
    % A lot of the parameters are shared among trials and conditions, so
    % for tidyness and for parfor keep them in a separate function.
    [cfgTEP,cfgTESS] = set_common_parameters();
    
    cfgTEP.toi            = [0,max(max(data.time{1}))]; % time of interest, make sure trials have the same length.
    
    % optimizing embedding
    cfgTEP.ragdim         = data.dim(1,1);       % criterion dimension
    cfgTEP.repPred        = round(size(data.trial{1},2)*(1/2));

    % results file name
    cfgTESS.fileidout  = fullfile(OutputDataPath,[DATA{f}.pp{1} DATA{1}.conditions(1,3)]);
    cfgTEP.outputpath = OutputDataPath;
    
    % How many iterations of boostrapping and which block-shuffling method?
    shuffle_n = 1e2; % 1e3

    % I'm not sure which of these is best. You might have to try all.
    cfgTESS.surrogatetype  = 'blockreverse4';
    %cfgTESS.surrogatetype  = 'dd_blockreverse1';
    %cfgTESS.surrogatetype  = 'dd_blockreverse2';

    %% Start! (This can literally take hours!)
    TE_rez_temp = dd_TEshufftest(cfgTEP,cfgTESS,data,strcat(cfgTEP.outputpath,DATA{f}.pp{1}),shuffle_n);
    TE_rez{f} = TE_rez_temp;
    
    for r = 1:2
        for n = 1:size(TE_rez_temp.TEmat_shuffle,2)
            TE_rez{f}.surr_test(r,n)=sum(squeeze(TE_rez_temp.TEmat_shuffle(r,n,:))<TE_rez_temp.TEmat(r,n))/size(TE_rez_temp.TEmat_shuffle,3)*100;
        end
    end
    % For quick check look at the TE and the surrogate TE distributions in the first few trials:
    % for r=1:2;for tr=1:8;subplot(2,8,(r-1)*8+tr);hist(squeeze(TE_rez{f}.TEmat_shuffle(r,tr,:)),1e1);hold on;plot(TE_rez{f}.TEmat(r,tr),0,'^r');xlim([-.05 .1]);hold off;end;end
end
% save(fullfile(OutputDataPath,['TE_rez_',datestr(date,'yyyy-mm-dd'),'.mat']),'TE_rez')
save(fullfile(OutputDataPath,['mid_processing_' datestr(date,'YYYYMMDD-hhmmss') '.mat']),'-v7.3')

return


%% To get a feel of the TE wrt to task performance.
% !! This is just an example! You need to extract the TE surrogate test 
% results in a way that aggregates the trials from the same participant in 
% the right order. Be careful! Right now the cells in DATA are broken both
% by participant and stimulus type.
for f=1:numel(DATA)
    figure(1)
    s=dlmread(fullfile(InputDataPath,DATA{f}.pp{1},'scores'),',',1,1);
    subplot(1,2,1)
    plot(TE_rez{f}.surr_test','-o')
    hold on
    plot(s(:,1)*1e2,'-s','linewidth',2)
    legend('TE_{Stim->User}','TE_{User->Stim}','Performance score')
    hold off
    subplot(1,2,2)
    scatter(s(:,1),TE_rez{f}.surr_test(1,:)','b')
    hold on
    scatter(s(:,1),TE_rez{f}.surr_test(2,:)','r')
    hold off
    fprintf('%3.0f\n',DATA{f}.conditions(1,3))
    pause
end

Surr = zeros(40,2,3);
Surrcount= zeros(40,1,3);
for f=1:numel(DATA)
    switch DATA{f}.conditions(1,3)
        case 10
            task = 1;
        case 25
            task = 2;
        case 30
            task = 3;
    end
    Surrcount(1:numel(TE_rez{f}.surr_test(1,:)'),1,task) = Surrcount(1:numel(TE_rez{f}.surr_test(1,:)'),1,task)+1;
    Surr(1:numel(TE_rez{f}.surr_test(1,:)'),1,task) = Surr(1:numel(TE_rez{f}.surr_test(1,:)'),1,task)+TE_rez{f}.surr_test(1,:)';
    Surr(1:numel(TE_rez{f}.surr_test(1,:)'),2,task) = Surr(1:numel(TE_rez{f}.surr_test(1,:)'),2,task)+TE_rez{f}.surr_test(2,:)';
end
for task=1:3
    Surr(:,1,task) = Surr(:,1,task)./Surrcount(:,1,task);
    Surr(:,2,task) = Surr(:,2,task)./Surrcount(:,1,task);
    Surr(:,3,task) = task;
end
Surr = reshape(permute(Surr,[1 3 2]),[],3);
for task=1:3
    subplot(1,3,task)
    plot(Surr(Surr(:,3)==task,1),'-ob')
    hold on
    plot(Surr(Surr(:,3)==task,2),'-or')
    hold off
end


for f=1:numel(DATA)
    % Inspect the raw trajectory and the TE outcomes in both directions.
    % For example, if the participant's amplitude is with a consistently 
    % larger amplitude than the stimulus (the participant is "driving" the
    % stimulus) you expect to see much larger values in the user->stim
    % surrogate test than in the stim->user.
    figure(2)
    for tr = 1:size(DATA{f}.trial,2)
        plot(DATA{f}.time{tr},DATA{f}.trial{tr}')
        fprintf('%6.2f,%6.2f,%6.2f\n',s(tr),TE_rez{f}.surr_test(:,tr)')
        pause
    end
end


%% Stats, trends with practice trials, etc.
% You need to pull the surrogate test scores from the structure and put them
% in an array that also has participant number, condition, trial number, etc.. 
% Such a table in the long format is to be used for stats and regressions.

