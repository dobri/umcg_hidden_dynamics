function tau = tau_quarter_period(x,t,varargin)

if numel(varargin)>0
    plotting_flag = varargin{1};
else
    plotting_flag = 0;
end

dt=mean(diff(t));
x=x-mean(x);
x=detrend(x);
x=x./max(x);
[~, pklocs]=findpeaks(x,'MinPeakHeight',.2,'MinPeakDistance',.5*1/dt);

if plotting_flag
    plot(t,x);hold on;plot(t(pklocs),x(pklocs),'+r');hold off;
end

tau=round(mean(diff(pklocs))/4);
end