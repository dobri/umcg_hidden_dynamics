% cd('C:\Users\ddotov\biomech\projects\umcg\data\wii_study\data')
cd('/home/ddotov/biomech/projects/umcg/data/wii_study/data/')
clear
pps = dir('PP*')';
A = table('Size',[3000,12],'VariableTypes',["int8","int8","string",'int8','string','logical','double','int32','double','double','double','double']);
A.Properties.VariableNames = {'PP','Day','Date','TaskCode','TrainingPhase','Visual','TrialDuration','SamplesNum','Score','C','tau','PitchError'};
c = 0;
for p = pps
    dag = [];
    dag{1} = dir(fullfile(p.folder,p.name,'Dag1','trial_log*'))';
    dag{2} = dir(fullfile(p.folder,p.name,'Dag2','trial_log*'))';
    % if isempty(dag{1});continue;end
    % if isempty(dag{2});continue;end
    for d = 1:2
        trial_counter = 0;
        S = dir(fullfile(dag{d}(1).folder,'scores*'));
        % corrupt scores file in pp 29
        if isempty(S);continue;end
        scores = readtable(fullfile(S(end).folder,S(end).name));
        %for f = dag{d}
        for s = 1:size(scores,1)
            fname = scores.rawDataFile{s};
            x = dlmread(fullfile(dag{d}(1).folder,fname),',',1,0);
            % There were corrupted raw data files in PP29, and empty
            % folders in PP 31 and 32. These participants were eliminated.
            % try x = dlmread(fullfile(dag{d}(1).folder,fname),',',1,0);catch me;disp([p d s]);keyboard;continue;end
            % x = readtable(fullfile(f.folder,fname));
            trial_counter = trial_counter + 1;
            c = c+1;
            A{c,1} = str2double(regexp(fname,'((?<=pp).*(?=_d))','match','once'));
            A{c,2} = str2double(regexp(fname,'(?<=_d)([^_]+)','match','once'));
            A{c,3} = string(regexp(fname,'((?<=log-).*(?=_task))','match','once'));
            A{c,4} = str2double(regexp(fname,'((?<=task).*(?=_aud))','match','once'));
            if d == 1
                if any(regexp(fname,'test'))
                    if trial_counter < 30
                        A{c,5} = "PreTest__";
                    else
                        A{c,5} = "PostTest_";
                    end
                else
                    A{c,5} = "Training_";
                end
            else
                A{c,5} = "Retention";
            end
            A{c,6} = str2double(regexp(fname,'((?<=vis).*(?=_eps))','match','once'))==1;
            A{c,7} = x(end,1);
            A{c,8} = size(x,1);
            A{c,9} = scores.score(s);
            A{c,10} = scores.cmax(s);
            A{c,11} = scores.tau(s);
            A{c,12} = scores.rmse(s);
        end
    end
end
A(c+1:end,:) = [];
save(['Scores_' char(datetime('now','Format','yyyy-MM-dd')) '.mat'],'A')
