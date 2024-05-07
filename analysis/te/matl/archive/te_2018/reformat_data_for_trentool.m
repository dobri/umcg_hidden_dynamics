function fileCell = reformat_data_for_trentool(raw_data_file,InputDataPath,variable,selected_lv)
fileCell=cell(1,32); % 8 trials x 4 blocks
total_duration_to_analyze = 160;
switch variable
    case 'lvpos'
        load(raw_data_file,'lvpos','lvtime','lvmeta','censure')
        raw_data=lvpos;
    case 'lvenergy'
        load(raw_data_file,'lvenergy','lvtime','lvmeta','censure')
        raw_data=lvenergy;
    case 'vel' % !
end

c=0;
for trial=1:8
    raw_data{trial}(:,:,censure)=[];
    for block=1:4
        c=c+1;
        in1=floor((block-1)*(total_duration_to_analyze/4*lvmeta.fr))+1;
        in2=floor(block*(total_duration_to_analyze/4*lvmeta.fr));
        data.trial{1}=squeeze(raw_data{trial}(in1:in2,selected_lv,:))';
        data.time{1}=lvtime{trial}(in1:in2)'-lvtime{trial}(in1)+1/lvmeta.fr;
        data.fsample=lvmeta.fr;
        for pp=1:size(raw_data{trial},3)
            data.label{1,pp}=['P' num2str(pp+(pp>=censure),'%02.f')];
        end
        fileCell{c} = ['RawData_LvPos' num2str(selected_lv) '_Trial' num2str(trial,'%02.f') '_Block' num2str(block,'%02.f') '.mat'];
        fileoutsave=fullfile(InputDataPath,fileCell{c});
        save(fileoutsave,'data')
    end
end
save(fullfile(InputDataPath,'fileCell_array.mat'),'fileCell')
