function DATA = convert_csv_to_ft_and_te(InputDataPath)
%%
%{
Loop through participant folders and trials to then process in batches.
Group the trials in a structure with levels for participant. Some analysis 
parameters will be constrained to be the same for all trials in the current
batch. Also useful for parallel processing which is needed because the 
analysis takes a lot of time.
+ optional:
Instead of wasting computation time and allowing for error, we can set
the embedding dimension based on our prior knowledge of the task dynamic.
m=2 for pure oscillaton (phase oscillator) and m=3 for nonlinear, chaotic.
The raw data is also pre-analyzed to get tau (for the embedding delay).
%}
ds = 2;

%% Prepare the data, arrange in the right format, get the conditions for the file name, etc..

% Assuming that each particpant's data is in a separate folder that is
% labeled in a consistent manner. Here for demonstration I just used data
% where I was the pilot and a copy of it (for testing parfor and
% consistency of surrogate testing).
folderCell = dir(fullfile(InputDataPath,'pp*'));

DATA = cell(1); % The number of cells in DATA will grow.
counter = 0;
for f=1:length(folderCell)
    S = readtable(fullfile(folderCell(f).folder,folderCell(f).name,'scores'));
    f_trials = dir(fullfile(folderCell(f).folder,folderCell(f).name,'trial*'));
    
    %     % How many different stimulus types in the given folder?
    %     tasks_vec = [];
    %     for tr = 1:numel(f_trials)
    %         tasks_vec(tr) = str2double(f_trials(tr).name(29:30));
    %     end
    %     tasks = unique(tasks_vec);
    %
    %     % Gather by stimulus type and then trials that are of the given type.
    %     for ta = 1:numel(tasks) % That's not needed here. Tasks are fully b/w.
    counter = counter + 1;
    counter_d = 0;
    for tr = 1:numel(f_trials)
        %             if tasks_vec(tr) == tasks(ta)
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
        
        % Check the printed output here to make sure that there aren't 
        % unusually short trials and other irregularities.
        fprintf('%8.2f,',cond_vec,max(t),1/mean(diff(t)));fprintf('\n')
        if max(t)<40;keyboard;end
        
        
        %% A little bit of pre-processing: smooth and detrend.
        %plot(x_stim,'-b');hold on
        %plot(x_user,'-m')
        x_stim = smooth(x_stim,10,'sgolay');
        x_user = smooth(x_user,10,'sgolay');
        % linear detrend the participant movement cause this can
        % mess up several of the parameters (inf corr time).
        b = [ones(size(x_user)) t]\x_user;
        x_user = x_user - b(2)*t;
        
        %plot(x_stim,'--b')
        %plot(x_user,'--m')
        %hold off;pause

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
        
        t=t(1:ds:end);
        x_stim=x_stim(1:ds:end);
        x_user=x_user(1:ds:end);
        
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
        %end
    end
    
    %% Just in case crop trials to same length. That's how fieldtrip wants them.
    maxl = Inf;
    for tr = 1:size(DATA{counter}.time,2)
        maxl = min([maxl size(DATA{counter}.time{tr},2)]);
    end
    for tr = 1:size(DATA{counter}.time,2)
        DATA{counter}.trial{tr} = DATA{counter}.trial{tr}(:, 1:maxl);
        DATA{counter}.time{tr} = DATA{counter}.time{tr}(:, 1:maxl);
    end
    %end
    fprintf('\n\n\n')
    %pause
    DATA{counter}.fsample = round(mean(DATA{counter}.fsample(counter_d)));
end

return