%addpath(('F:\nonlinear_interdependence\'))
flbase='~/c3/experiment1/data';
raw_output_folder='~/Desktop';
cd(flbase)
folders=dir('pp*');
S=[];
for k=1:length(folders)
    fprintf('%s\n',folders(k).name);
    s=loop_analysis_te(flbase,[folders(k).name '/'],raw_output_folder);
    % We need the tau and dim.
    %figure_nlinter_and_cc_rmse([flbase folders(k).name '\'])
    S=[S;s];
end
