function [Index_Conditions, Condition_tags] = Get_Diagnosis_from_CSV(results_dir)

load(fullfile(results_dir,'LEIDA_Eigenvectors.mat'),'Data_info')

CSV_table=importdata('DONAU_ADNI_fMRI_351subs.csv');
CSV_table=CSV_table.textdata;

Index_Conditions=zeros(1,size(Data_info,1));

for s=1:size(Data_info,1)
    sub_name=Data_info(s).name;
    sub_name=sub_name(5:12);
    sub_name=[sub_name(1:3) '_S_' sub_name(5:8)];
    Data_info(s).subject=sub_name;
    idx = find(contains(CSV_table(:,1), sub_name));
    Data_info(s).diagnose=char(CSV_table(idx,17));
    Data_info(s).diagnoseIdx=str2double(char(CSV_table(idx,16)));
    Index_Conditions(s)=str2double(char(CSV_table(idx,16)));
end

for cond=1:max(Index_Conditions)
    s=Index_Conditions==cond;
    Condition_tags{cond}=Data_info(s).diagnose;
end


save LEIDA_Eigenvectors.mat Data_info -append
