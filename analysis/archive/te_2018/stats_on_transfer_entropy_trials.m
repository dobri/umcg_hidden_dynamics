load ~/c3/experiment1/transfer_entropy/out/TGA_results_trials_and_data_and_params2017-10-01.mat

S=zeros(48,4,2);
x=csvread('~/c3/experiment1/analyzed/Deltas-02-Jun-2017.csv',1,0);
group=x(1:3:end,2);
S(:,1,1)=group;
S(:,1,2)=group;

Sdelta=zeros(48,4,1);
Sdelta(:,1)=group;


Strial=cell(1,3);
for test=1:3
    Strial{test}=cell(1,2);
end

for pp=1:48
    for test=1:3
        trial=(pp-1)*3+test;
        %disp(trial)
        surtest=[];
        for direction=1:2
            for tt=1:size(TGA_results_trials{trial}.TEmat_shuffle,2)
                % percentile
                surtest(tt,direction)=...
                    sum(squeeze(TGA_results_trials{trial}.TEmat_shuffle(direction,tt,:))...
                    <TGA_results_trials{trial}.TEmat(direction,tt))...
                    ./numel(squeeze(TGA_results_trials{trial}.TEmat_shuffle(direction,tt,:)));
            end
        end
        S(pp,test+1,1)=mean(surtest(:,1));
        S(pp,test+1,2)=mean(surtest(:,2));
        if test==3
            gr = S(trial/3,1,1);
            temp=[surtest;nan(50-length(surtest),2)];
            for direction=1:2
                Strial{gr}{direction}=horzcat(Strial{gr}{direction},temp(:,direction));
            end
        end
        if test==2
            gr = S((trial+1)/3,1,1);
            disp([pp gr trial])
            if mod(size(surtest,1),2)==1
                delta1=surtest(1:3,1);
                delta2=surtest(4:7,1);
            elseif size(surtest,1)==4
                delta1=surtest(1:2,1);
                delta2=surtest(3:4,1);
            else
                delta1=surtest(1:size(surtest,1)/2,1);
                delta2=surtest((size(surtest,1)/2+1):end,1);
            end
            Sdelta(pp,3)=mean(delta1);
            Sdelta(pp,4)=mean(delta2);
            Sdelta(pp,2)=(Sdelta(pp,4)-Sdelta(pp,3))./Sdelta(pp,3)*100;
        end
        %         plot(surtest)
        %         pause
    end
end



%{

[~,p,~,stats]=ttest(Sdelta(Sdelta(:,1)==1,2),0,.05,'right')
p =
    0.0961
stats = 
  struct with fields:

    tstat: 1.3654
       df: 15
       sd: 35.7630

[~,p,~,stats]=ttest(Sdelta(Sdelta(:,1)==2,2),0,.05,'right')
p =
    0.0809
stats = 
  struct with fields:

    tstat: 1.4773
       df: 14
       sd: 39.5001

[~,p,~,stats]=ttest(Sdelta(Sdelta(:,1)==3,2),0,.05,'right')
p =
    0.0248
stats = 
  struct with fields:

    tstat: 2.1240
       df: 16
       sd: 13.9117

%}