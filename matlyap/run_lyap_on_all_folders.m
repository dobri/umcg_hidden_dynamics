clear
addpath('/home/dobri/logos/c3/umcg_hidden_dynamics/matlyap')
flbase = '/media/dobri/disk2/c3/umcg/data/training_processed';
cd(flbase)
folders=dir('pp*');
S=[];
if isempty(dir('scores_lyapunov*'))
    labels = [{'pp'},'trial','condition','period',{'LyapsTutor'},'LyaplTutor','LyapsTrainee','LyaplTrainee'];
    fid=fopen(fullfile(flbase,['scores_lyapunov-',datestr(date),'.csv']),'a');
    for kk=1:length(labels);fprintf(fid,'%s;',labels{kk});end;fprintf(fid,'\n');
    fclose(fid);
end
for k = 1:length(folders)
    fprintf('%s\n',folders(k).name);
    s=loop_analysis_lyap_divergence(flbase,[folders(k).name '/'],0,0);
    S=[S;s];
    if 1
        fid=fopen(fullfile(flbase,['scores_lyapunov-',datestr(date),'.csv']),'a');
        for kk=1:size(s,1);fprintf(fid,'%.4f;',s(kk,:));fprintf(fid,'\n');end
        fclose(fid);
    end
end