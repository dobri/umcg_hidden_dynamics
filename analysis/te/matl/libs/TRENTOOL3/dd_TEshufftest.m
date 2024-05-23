function TE_shufftest = dd_TEshufftest(cfgTEP,cfgTESS,data,fileidout,shuffle_n)
% This function has to be supplied with the folowing INPUT:
%
% data      -    the data in the fieldtrip raw data format; see TEprepare for more
%           details
% cfgTEP    -    a configuration structure that has all the fields required for
%           TEprepare; see TEprepare for more details. The only difference is that instead of
%           cfg.predictime_u which is a single number, cfg.predicttimemin_u
%           cfg.predicttimemax_u,cfg.predicttimestepsize
%           have to be supplied, indicating the minimum and maximum prediction time of
%           inteterest and the stepsize of the resolution
% cfgTSS    -   a configuration structure that has all the fields required
%           for TEsurrogatestats
%
% The OUTPUT TE_shufftest is a structure containing estimated TE values and 
% results for surrogate testing
% using the optimal interaction delay for each channel combination?

t_total = tic;

%% check input
% cfgTEP.toi=toi;
% cfgTEP.repPred=repPred;
% cfgTEP.channel=channel;
cfgTESS.fileidout=fileidout;

%% define logging levels
LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;

if ~isfield(cfgTEP, 'verbosity'), cfgTEP.verbosity = 'info_minor'; end;

%% checks and parameter preparations

% cfg.predicttimevec_u supplied ?
if ~isfield(cfgTEP,'predicttimemax_u')
    fprintf('\n')
    error('TRENTOOL ERROR: No cfgTEP.predicttimemax_u specified - see HELP InteractionDelayReconstruction_calculate for more information');
end
if ~isfield(cfgTEP,'predicttimemin_u')
    fprintf('\n')
    error('TRENTOOL ERROR: No cfgTEP.predicttimemin_u specified - see HELP InteractionDelayReconstruction_calculate for more information');
end
if ~isfield(cfgTEP,'predicttimestepsize')
    fprintf('\n')
    error('TRENTOOL ERROR: No cfgTEP.predicttimestepsize specified - see HELP InteractionDelayReconstruction_calculate for more information');
end

predicttimevec_u=cfgTEP.predicttimemin_u:cfgTEP.predicttimestepsize:cfgTEP.predicttimemax_u;

if ~(predicttimevec_u(end)==cfgTEP.predicttimemax_u)
    predicttimevec_u(end+1)=cfgTEP.predicttimemax_u; % make sure the last intended ptredictiontime is also investigated
end
% all other checks are left to the subsidiary functions


%% TEprepare part
msg = '################### PREPARING DATA FOR TE ANALYSIS';
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MAJOR);

cfgTEP.predicttime_u = 0;%cfgTEP.predicttimemax_u; %cfgTEP.predicttimemax_u;  % fix config
dataprep = TEprepare(cfgTEP,data);
clear data;


%% Find optimal interaction delays. Or not.
msg = '################### OPTIMIZING INFORMATION TRANSFER DELAY';
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MAJOR);

if cfgTEP.predicttimemin_u == cfgTEP.predicttimemax_u
    dataprep.TEprepare.u_in_ms = repmat(cfgTEP.predicttime_u,size(cfgTEP.channel,2)^2-size(cfgTEP.channel,2),1);
    dataprep.TEprepare.u_in_samples = repmat(round(cfgTEP.predicttime_u/1000*dataprep.fsample),size(cfgTEP.channel,2)^2-size(cfgTEP.channel,2),1);
    dataprep.TEprepare.cfg.predicttime_u = cfgTEP.predicttime_u;
else
    dataprep = TEfindDelay(predicttimevec_u,cfgTESS,dataprep);
end
cfgTESS.embedsource = 'yes';

cfgTESS.extracond = 'none';


%% calulate statistics with optimal u for individual channels
msg = '################### ESTIMATING TRANSFER ENTROPY WITH OPTIMIZED PARAMETERS';
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MAJOR);

cfgTESS.fileidout=strcat(cfgTESS.fileidout,'_RAG4_TGA_opt_u');

TE_shufftest=dd_TEsurrogatestats(cfgTESS,dataprep,shuffle_n);

%%
t=toc(t_total);
msg = sprintf( ...
    'Thank you for using this transfer entropy tool!\n\nTRANSFER ENTROPY CALCULATION ENDED: %s \nCALCULATION TOOK %.0f MINUTES (%.0f SECONDS)', ...
    datestr(now), t/60, t);
TEconsoleoutput(cfgTEP.verbosity, msg, LOG_INFO_MAJOR);