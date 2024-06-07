clear

%% set software paths
addpath('~/logos/umcg_hidden_dynamics/analysis/maxlyap/matl/')

%% define data paths
base_folder = '~/biomech/projects/side_projects/umcg/data/';
% base_folder = 'C:\Users\ddotov\biomech\projects\side_projects\umcg\data\';
OutputDataPath = fullfile(base_folder,'results/');
InputDataPath  = fullfile(base_folder,'raw_data/');
plotting = 0;

cd(fullfile(base_folder))
data_aggregated_already = 1;
if data_aggregated_already
    load(fullfile(InputDataPath,'DATA_2024-06-06.mat'))
else
    DATA = load_raw_in_DATA_per_task(InputDataPath);
end


%%
T = table('Size',[3000,16],'VariableTypes',{'string','int8','int8','string',...
    'int8','string','string','int8',...
    'double','double','double','double','double','double','double','double'});
T.Properties.VariableNames = {'FileName','PP','Day','Date',...
    'TaskCode','TrainingPhase','TrainingPhaseLabel','Visual',...
    'Score','C','tau','PitchError',...
    'MaxLyapStimShort','MaxLyapStimLong','MaxLyapUserShort','MaxLyapUserLong'};
counter = 0;
for task = 1:numel(DATA)
    data = DATA{task};
    for tr = 1:size(data.trial,2)
        x_user = data.trial{tr}(2,:)';
        dt = 1/data.fsample;
        meanperiod = data.tau(1,tr)*4;
        if DATA{task}.task{1,tr} < 20
            x = data.generator{tr}';
            x = smooth(x,5,'lowess');
            v = [diff(x)./dt; 0];
            v = smooth(v,5,'lowess');
            acc = [diff(v)./dt; 0];
            acc = smooth(acc,5,'lowess');
            s = [x v acc];
        else
            s = data.generator{1}';
        end
        [lyaptutors,lyaptutorl] = max_lyapunov_script(s,meanperiod,dt,plotting,[]);
        [lyaplearners,lyaplearnerl] = max_lyapunov_script(x_user,meanperiod,dt,plotting,[]);

        fprintf('%s\n',DATA{task}.trial_name{1,tr})
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
        T{counter,13:16} = [lyaptutors lyaptutorl lyaplearners lyaplearnerl];
    end
end
T(counter+1:end,:) = [];


writetable(T,fullfile(OutputDataPath,['Scores_MaxLyap_' char(datetime('now','Format','yyyy-MM-dd')) '.csv']))
%{
%}