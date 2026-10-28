function LEiDA_Design(P)
%
% Second step of LEiDA (called from LEiDA_Main): for the study design
% P.Design, assign each scan to a condition, compute the fractional
% occupancy of each FC state for all K, compare it between conditions with
% permutation tests, and plot the results: p-values and effect sizes across
% K, barplots and the centroid pyramid. All results are saved in the
% folder of the design (P.design_dir).
%
% INPUT:
% P         parameters set in LEiDA_Main (after LEiDA_Start)
%
% Tutorial: README.md
% Version:  V2.0, October 2026
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

D = P.Design;
if ~exist(P.design_dir,'dir')
    mkdir(P.design_dir);
end
disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STUDY DESIGN: ' P.Design_name ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
disp(' ')

if strcmpi(D.Condition_source,'none')
    Condition_tags = {};
else
    % Assign each scan to a condition
    if ~isfield(D,'Condition_column')
        D.Condition_column = '';
    end
    [Index_Conditions,Condition_tags] = LEiDA_conditions(P.leida_results,D.Condition_source,...
        D.Condition_tags,D.Metadata_file,D.Metadata_subject_column,D.Metadata_session_column,...
        D.Condition_column,P.design_dir);
    if isfield(D,'Condition_names') && ~isempty(D.Condition_names)
        if numel(D.Condition_names) ~= numel(Condition_tags)
            error('LEiDA:conditions','Condition_names must have one name per condition in Condition_tags');
        end
        Condition_tags = D.Condition_names(:)';
    end

    % Compute the fractional occupancy and perform hypothesis tests
    LEiDA_stats_FracOccup(P.leida_results,Index_Conditions,Condition_tags,D.Paired_tests,...
                          D.n_permutations,D.n_bootstraps,D.Stats_level,P.design_dir);

    % Print the 10 largest effect sizes of each pair of conditions, with
    % their 95% bootstrap confidence intervals and p-values
    S = load([P.design_dir 'LEiDA_Stats_FracOccup.mat'],'effectsize_signed','effectsize_ci','P_pval2sided','rangeK');
    [cond2,cond1] = find(triu(ones(numel(Condition_tags)),1)');
    for pair = 1:numel(cond1)
        g_pair = reshape(S.effectsize_signed(pair,:,:),[],1);
        [~,idx] = sort(abs(g_pair),'descend');
        idx = idx(1:min(10,numel(idx)));
        [k,c] = ind2sub([numel(S.rangeK) S.rangeK(end)],idx);
        ci = reshape(S.effectsize_ci(pair,:,:,:),[],2);
        p = reshape(S.P_pval2sided(pair,:,:),[],1);
        disp(['Largest effect sizes, ' Condition_tags{cond1(pair)} ' vs ' Condition_tags{cond2(pair)} ...
              ' (g > 0: ' Condition_tags{cond1(pair)} ' higher):'])
        disp(table(S.rangeK(k)',c,g_pair(idx),ci(idx,1),ci(idx,2),p(idx),...
                   'VariableNames',{'K','c','g','CI_low','CI_high','p'}))
    end

    % Generate and save the p-value and barplot plots for fractional occupancy
    Plot_FracOccup(P.design_dir)
end

% Plot the centroids obtained using LEiDA, their overlap with Yeo nets and
% their fractional occupancy in each condition
Plot_Centroid_Pyramid(P.leida_results,Condition_tags,P.Parcellation,P.N_areas,P.CortexDirection,[],P.design_dir)

% Save the parameters of this design for reference
save([P.design_dir 'LEiDA_Parameters.mat'],'-struct','P');
