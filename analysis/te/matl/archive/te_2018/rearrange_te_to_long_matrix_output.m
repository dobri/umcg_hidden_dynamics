load TGA_results_trials_and_data_and_params2017-10-01.mat

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
                surtest(tt,direction)=...
                    sum(squeeze(TGA_results_trials{trial}.TEmat_shuffle(direction,tt,:))...
                    <TGA_results_trials{trial}.TEmat(direction,tt))...
                    ./numel(squeeze(TGA_results_trials{trial}.TEmat_shuffle(direction,tt,:)));
                SURTEST{pp,test}=surtest;
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
            Delta=(mean(delta2)-mean(delta1))./mean(delta1)*100;
            Sdelta(pp,2)=Delta;
        end
        %plot(surtest)
        %pause
    end
end

return

xx=csvread('~/c3/experiment1/analyzed/master_data_file.csv',1,0);

% Redistribute back to the simple long format of the master output csv file
X=[];
for pp=1:48
    pp
    for type = 1:5
        type
        switch type
            case 1
                index=logical((xx(:,1)==pp).*(xx(:,12)==1).*(xx(:,14)<=2));
                X=vertcat(X,[xx(index,1:14),SURTEST{pp,1}(1:sum(index),:)]);
            case 2
                index=logical((xx(:,1)==pp).*(xx(:,12)==1).*(xx(:,14)==3));
                X=vertcat(X,[xx(index,1:14),SURTEST{pp,2}(1:sum(index),:)]);
            case 3
                index=logical((xx(:,1)==pp).*(xx(:,12)==0).*(xx(:,13)==0));
                X=vertcat(X,[xx(index,1:14),SURTEST{pp,3}(:,:)]);
            case 4
                index=logical((xx(:,1)==pp).*(xx(:,13)==1).*(xx(:,14)<=2));
                X=vertcat(X,[xx(index,1:14),SURTEST{pp,1}((end-sum(index)+1):end,:)]);
            case 5
                index=logical((xx(:,1)==pp).*(xx(:,13)==1).*(xx(:,14)==3));
                X=vertcat(X,[xx(index,1:14),SURTEST{pp,2}((end-sum(index)+1):end,:)]);
        end
    end
end

% fid=fopen('temp.csv','a');for r=1:size(X,1);fprintf(fid,'%10.4f,',X(r,:));fprintf(fid,'\n');end;fclose(fid);