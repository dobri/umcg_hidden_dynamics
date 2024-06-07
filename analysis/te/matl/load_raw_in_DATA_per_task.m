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

% Prepare the data, arrange in the right format, get the conditions for the file name, etc..
cd('~/biomech/projects/side_projects/umcg/data/raw_data/')
clear
pps = dir('PP*')';
A = table('Size',[3000,14],'VariableTypes',{'string','int8','int8','string','int8','string','int8','double','int32','double','double','double','double','int8'});
A.Properties.VariableNames = {'FileName','PP','Day','Date','TaskCode','TrainingPhase','Visual','TrialDuration','SamplesNum','Score','C','tau','PitchError','TrainingPhaseNumber'};

tasks = [11 12 20 21 10 25 30];
DATA = cell(1,numel(tasks));

counter_trial = zeros(1,7);
c = 0;
ppn = 0;
for p = pps
    ppn = ppn + 1;
    dag = [];
    dag{1} = dir(fullfile(p.folder,p.name,'Dag1','trial_log*'))';
    dag{2} = dir(fullfile(p.folder,p.name,'Dag2','trial_log*'))';
    for d = 1:2
        trial_counter = 0;
        S = dir(fullfile(dag{d}(1).folder,'scores*'));
        if isempty(S);continue;end % corrupt scores file
        scores = readtable(fullfile(S(end).folder,S(end).name));
        for s = 1:size(scores,1)
            fname = scores.rawDataFile{s};
            x = dlmread(fullfile(dag{d}(1).folder,fname),',',1,0);
            t = x(:,1); % time
            if(x(end,1)<59);continue;end

            % Trial conditions and parameters. Same as in performance loop
            trial_counter = trial_counter + 1;
            c = c+1;
            A{c,1} = string(fname);
            A{c,2} = str2double(regexp(fname,'((?<=pp).*(?=_d))','match','once'));
            A{c,3} = str2double(regexp(fname,'(?<=_d)([^_]+)','match','once'));
            A{c,4} = string(regexp(fname,'((?<=log-).*(?=_task))','match','once'));
            A{c,5} = str2double(regexp(fname,'((?<=task).*(?=_aud))','match','once'));
            if d == 1
                if any(regexp(fname,'test'))
                    if trial_counter < 30
                        A{c,6} = "PreTest__";
                        A{c,14} = -1;
                    else
                        A{c,6} = "PostTest_";
                        A{c,14} = 1;
                    end
                else
                    A{c,6} = "Training_";
                    A{c,14} = 0;
                end
            else
                A{c,6} = "Retention";
                A{c,14} = 2;
            end
            A{c,7} = str2double(regexp(fname,'((?<=vis).*(?=_eps))','match','once'))==1;
            A{c,5} = A{c,5} + A{c,7};
            A{c,8} = x(end,1);
            A{c,9} = size(x,1);
            A{c,10} = scores.score(s);
            A{c,11} = scores.cmax(s);
            A{c,12} = scores.tau(s);
            A{c,13} = scores.rmse(s);

            % Don't search but decide the embedding dimension: 2 or 3,
            % depending on the task space (oscillator or nonlinear sys)
            if A{c,5} >= 20
                m_stim = 3;
                m_pp = 3;
            else
                m_stim = 3.0;
                m_pp = 3.0;
            end

            % Prepare for TE
            if A{c,7} == 1 % visual modality condition.
                x_stim = x(:,16); % the stimulus
                x_user = x(:,17); % the participant
            else
                x_stim = x(:,11); % the stimulus
                x_user = x(:,12); % the participant
            end
            if A{c,5} >= 20
                x_gen = x(:,7:9); % the generator
            else
                x_gen = x(:,8);
            end

            % A little bit of pre-processing: smooth and detrend.
            x_stim = smooth(x_stim,10,'sgolay');
            x_user = smooth(x_user,10,'sgolay');
            for kk = 1:size(x_gen,2)
                x_gen(:,kk) = smooth(x_gen(:,kk),10,'sgolay');
            end
            % linear detrend the participant movement because this can
            % mess up several of the parameters (inf corr time).
            b = [ones(size(x_user)) t]\x_user;
            x_user = x_user - b(2)*t;

            % Resample to an uniform timeframe.
            Fs = 50;
            [x_user, ~] = resample(x_user, t, Fs);
            [x_stim, ~] = resample(x_stim, t, Fs);
            [x_gen, t] = resample(x_gen, t, Fs);

            % Crop the first 10 seconds of the trial and the last 59+ seconds of the trial.
            x_user = x_user((t>10) & (t<=59),:);
            x_stim = x_stim((t>10) & (t<=59),:);
            x_gen = x_gen((t>10) & (t<=59),:);
            t = t((t>10) & (t<=59),:); % Crop the last 59+ seconds of the trial!

            if sum(isnan(x_stim)) > 0 ; keyboard ; end


            t = t - min(t); % for simplicity, shift back to zero.
            fprintf('%10.0f%10.4f%10.0f\n',[min(t) max(t) numel(t)])

            x_stim = (x_stim-mean(x_stim))./std(x_stim);
            x_user = (x_user-mean(x_stim))./std(x_stim);

            % Estimate tau even though for now we do not use it. :/
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

            celln = find(tasks==A{c,5});
            counter_trial(celln) = counter_trial(celln) + 1;
            trialn = counter_trial(celln);

            DATA{celln}.trial_name{trialn} = A{c,1};
            DATA{celln}.pp{trialn} = A{c,2};
            DATA{celln}.day{trialn} = A{c,3};
            DATA{celln}.date{trialn} = A{c,4};
            DATA{celln}.task{trialn} = A{c,5};
            DATA{celln}.training_phase{trialn} = A{c,14};
            DATA{celln}.training_phase_label{trialn} = A{c,6};
            DATA{celln}.visual{trialn} = A{c,7};
            DATA{celln}.duration{trialn} = A{c,8};

            DATA{celln}.scores(:,trialn) = A{c,10:13}';
            DATA{celln}.pp{trialn} = A{c,2};
            DATA{celln}.file_name{trialn} = A{c,1};
            DATA{celln}.conditions(:,trialn) = A{c,[5 7 14]}';
            DATA{celln}.time{trialn} = t';
            DATA{celln}.fsample = 1/mean(diff(t));
            DATA{celln}.label = {'Stim','User'};
            DATA{celln}.tau(:,trialn) = [tau_stim; tau_pp];
            DATA{celln}.dim(:,trialn) = [m_stim; m_pp];
            DATA{celln}.trial{trialn} = [x_stim x_user]';
            DATA{celln}.generator{trialn} = x_gen';
        end
    end
end
A(c+1:end,:) = [];

% save(['DATA_' char(datetime('now','Format','yyyy-MM-dd')) '.mat'],'DATA','A','-v7.3')
