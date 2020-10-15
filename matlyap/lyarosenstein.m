function d = lyarosenstein(x,m,tau,meanperiod,maxiter)
% d:divergence of nearest trajectoires
% x:signal
% tao:time delay
% m:embedding dimension
N=length(x);
M=N-(m-1)*tau;
Y=psr_deneme(x,m,tau);
neardis=nan(M,1);
nearpos=nan(M,1);
d=nan(maxiter,1);
for k=1:M
    x0=ones(M,1)*Y(k,:);
    distance=sqrt(sum((Y-x0).^2,2));
    for j=1:M
        if abs(j-k)<=meanperiod
            distance(j)=1e10;
        end
    end
    [neardis(k),nearpos(k)]=min(distance);
    %if mod(k,1000)==0;fprintf('%6.2f%s',k/M*100,'%');end
    %if mod(k,10000)==0;fprintf('\n');end
end
fprintf('\n');
for k=1:maxiter
    maxind=M-k;
    evolve=0;
    pnt=0;
    for j=1:M
        if j<=maxind && nearpos(j)<=maxind
            dist_k=sqrt(sum((Y(j+k,:)-Y(nearpos(j)+k,:)).^2,2));
            if dist_k~=0
                evolve=evolve+log(dist_k);
                pnt=pnt+1;
            end
        end
    end
    if pnt > 0
        d(k)=evolve/pnt;
    else
        d(k)=0;
    end
    if mod(k,100)==0;fprintf('%6.2f%s',k/maxiter*100,'%');end
    if mod(k,1000)==0;fprintf('\n');end
end
fprintf('\n')
%% LLE Calculation
%fs=128;%sampling frequency
%tlinear=15:78;
%F = polyfit(tlinear,d(tlinear),1);
%lle = F(1)*fs;
function Y=psr_deneme(x,m,tao,npoint)
%Phase space reconstruction
%x : time series
%m : embedding dimension
%tao : time delay
%npoint : total number of reconstructed vectors
%Y : M x m matrix
% author:"Merve Kizilkaya"
N=length(x);
if nargin == 4
    M=npoint;
else
    M=N-(m-1)*tao;
end
Y=zeros(M,m);
for i=1:m
    Y(:,i)=x((1:M)+(i-1)*tao)';
end