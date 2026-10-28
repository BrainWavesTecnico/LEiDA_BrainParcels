function Plot_FracOccup(data_dir)
%
% Plot the results from the hypothesis tests obtained from comparing the
% mean fractional occupancy between conditions
%
% INPUT:
% data_dir      directory with the results from running the hypothesis
%               tests on the mean fractional occupancy of FC states
%
% OUTPUT:
% Fig1          plot of the two-sided p-values obtained across K from the
%               comparison of the mean fractional occupancy of FC states
%               for each pair of conditions
% Fig2          barplot of the mean fractional occupancy of each FC state
%               for all conditions
% Fig3          barplot of the mean fractional occupancy of each FC state
%               for each pair of conditions
% Fig4          plot of the Hedge's effect size obtained across K from the
%               comparison of the fractional occupancy of FC states
%               for each pair of conditions
%
% Authors: Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with results for the fractional occupancy  (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';

% Load required data:
load([data_dir file_P],'cond', 'P', 'P_pval2sided', 'Index_Conditions', 'effectsize', 'rangeK', ...
     'Subjects', 'stats_level');
if ~exist('stats_level','var')
    stats_level = 'scan';
    Subjects = (1:numel(Index_Conditions))';
end

% Mean and SEM of the fractional occupancy in each condition, over the
% same observations used in the tests (subjects or scans)
[mean_P,sem_P] = FracOccup_mean_sem(P,Index_Conditions,Subjects,stats_level,numel(cond));

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Colors of the significance levels: p<0.05, p<0.05/K and p<0.05/sum(K)
darkRed = [0.6, 0, 0];
darkGreen = [0.3, 0.5, 0];
darkBlue = [0.07, 0.62, 1.00];

% Create a map to place the figures
pos = power(n_Cond-1,2);
index_fig = reshape(1:pos, n_Cond-1, n_Cond-1).';
subplot_map = ones(n_Cond-1);
subplot_map = triu(subplot_map).';

% Define indices for conditions in subplot
subplot_indices = find(subplot_map);

% Code below generates the figures relative to the fractional occupancy
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%% FRACTIONAL OCCUPANCY FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

% Possible pairs of conditions comparisons
condRow = zeros(1,n_Cond*(n_Cond-1)/2);
condCol = zeros(1,n_Cond*(n_Cond-1)/2);
cond_pair = 1;
disp('Pairs of conditions compared:')
for cond1 = 1:n_Cond-1
    for cond2 = cond1+1:n_Cond
        condRow(cond_pair) = cond1;
        condCol(cond_pair) = cond2;
        disp([num2str(cond{cond1}) ' - ' num2str(cond{cond2})])
        cond_pair = cond_pair + 1;
    end
end

% Plot of two-sided p-values from fractional occupancy hypothesis tests
Fig1 = figure('Position', [307   323   572   344]);
disp(' ')
disp('Plotting two-sided p-values from hypothesis (permutation) tests:')
for s_ind = 1:length(subplot_indices)
    
    % disp(['- ' cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}])
    
    subplot_ind = subplot_indices(s_ind);
    
    % Subplot adds plots following the rows (not the columns)
    subplot(size(subplot_map,1),size(subplot_map,2),index_fig(subplot_ind))
    
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'--','Color',darkRed,'LineWidth',1.5)
    hold on
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'--','Color',darkGreen,'LineWidth',1.5)
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'--','Color',darkBlue,'LineWidth',1.5)
    
    for k = 1:length(rangeK)
        for c = 1:rangeK(k)
            if P_pval2sided(s_ind,k,c) > 0.05
                semilogy(rangeK(k),P_pval2sided(s_ind,k,c),'.k','Markersize',8);
                % text(rangeK(k),P_pval2sided(s_ind,k,c),[' ' num2str(c)])
            end
            if P_pval2sided(s_ind,k,c) < 0.05 && P_pval2sided(s_ind,k,c) > (0.05/rangeK(k))
                semilogy(rangeK(k),P_pval2sided(s_ind,k,c),'+','Color',darkRed,'Markersize',8);
                % text(rangeK(k),P_pval2sided(s_ind,k,c),[' ' num2str(c)])
            end
            if P_pval2sided(s_ind,k,c) < (0.05/rangeK(k)) && P_pval2sided(s_ind,k,c) > (0.05/sum(rangeK))
                semilogy(rangeK(k),P_pval2sided(s_ind,k,c),'o','Color',darkGreen,'Markersize',8);
                % text(rangeK(k),P_pval2sided(s_ind,k,c),[' ' num2str(c)])
            end
            if P_pval2sided(s_ind,k,c) <= (0.05/sum(rangeK))
                semilogy(rangeK(k),P_pval2sided(s_ind,k,c),'*','Color',darkBlue,'Markersize',8);
                % text(rangeK(k),P_pval2sided(s_ind,k,c),[' ' num2str(c)])
            end
        end
    end
     
    title([cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}],'interpreter','none')
    ylabel([{'Two-sided {\itp}-value'},{'Fractional Occupancy'}])
    xlabel('Number of FC States K')
    xticks([2 4 6 8 10 12 14 16 18 20])
    xlim([rangeK(1)-1 rangeK(end)+1])
    box off
    ylim([10.^(-7.5) 1])
    set(gca,'YTick',10.^(-7:2:0))
end
Save_png(Fig1, fullfile(data_dir, 'FracOccup_pvalues.png'));
saveas(Fig1, fullfile(data_dir, 'FracOccup_pvalues.fig'),'fig');
disp('- Plot successfully saved as FracOccup_pvalues');

% Barplot of mean fractional occupancy for each condition across K (all),
% colored if there is a significant difference between any pair of
% conditions (smallest p-value over all pairs, Bonferroni-corrected for the
% number of pairs)
if n_Cond > 2
    n_compare = n_Cond*(n_Cond-1)/2;
    P_pval_any = reshape(min(min(P_pval2sided,[],1)*n_compare,1),size(P,2),size(P,3));
    Fig2 = figure('Position', get(0, 'Screensize'));
    disp(' ')
    disp('Plotting barplot of the mean fractional occupancy across conditions and K:')
    L = Plot_bar_pyramid(mean_P,sem_P,P_pval_any,rangeK,cond,1:n_Cond);
    annotation('textbox',[L.left 0 0.95 0.2*L.bottom],'FitBoxToText','off',...
        'VerticalAlignment','middle','EdgeColor','none','Interpreter','none','FontSize',10,...
        'String',['Fractional occupancy (mean ' char(177) ' SEM over ' stats_level 's). Color: difference ' ...
                  'between any pair of conditions (smallest p over ' num2str(n_compare) ' pairs x ' num2str(n_compare) ...
                  '): dark red p<0.05, dark green p<0.05/K, dark blue p<0.05/' num2str(sum(rangeK))])
    Save_png(Fig2, fullfile(data_dir, 'FracOccup_Barplot_Allconditions.png'));
    saveas(Fig2, fullfile(data_dir, 'FracOccup_Barplot_Allconditions.fig'),'fig');
    disp('- Plot successfully saved as FracOccup_Barplot_Allconditions');
end

% Barplot of mean fractional occupancy for each pair of conditions
disp(' ')
disp('Plotting barplots of mean fractional occupancy for each pair of conditions:')
n_compare = n_Cond*(n_Cond-1)/2;
for i = 1:n_compare
    Fig3 = figure('Position', get(0, 'Screensize'),'Name', ['Occupancy ' cond{condRow(i)} ' vs ' cond{condCol(i)}]);
    L = Plot_bar_pyramid(mean_P,sem_P,reshape(P_pval2sided(i,:,:),size(P,2),size(P,3)),rangeK,cond,[condRow(i) condCol(i)]);
    annotation('textbox',[L.left 0 0.95 0.2*L.bottom],'FitBoxToText','off',...
        'VerticalAlignment','middle','EdgeColor','none','Interpreter','none','FontSize',10,...
        'String',['Fractional occupancy (mean ' char(177) ' SEM over ' stats_level 's), ' cond{condRow(i)} ...
                  ' vs ' cond{condCol(i)} ': dark red p<0.05, dark green p<0.05/K, dark blue p<0.05/' num2str(sum(rangeK))])
    Save_png(Fig3, fullfile(data_dir, ['FracOccup_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)} '.png']));
    saveas(Fig3, fullfile(data_dir, ['FracOccup_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)} '.fig']),'fig');
    disp(['- Plot successfully saved as FracOccup_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)}]);
end

% Plot of Hedge's effect size from fractional occupancy hypothesis tests
disp(' ')
disp('Plotting Hedge''s effect size from hypothesis (permutation) tests:')
Fig4 = figure('Position', [307   323   572   344]);
for s_ind = 1:length(subplot_indices)
    
    subplot_ind = subplot_indices(s_ind);
    
    % Subplot adds plots following the rows (not the columns)
    subplot(size(subplot_map,1),size(subplot_map,2),index_fig(subplot_ind))
    
    plot([rangeK(1)-1 rangeK(end)+1],[0.8 0.8],'--','Color',darkBlue,'LineWidth',1.5);   
    hold on;
    plot([rangeK(1)-1 rangeK(end)+1],[0.5 0.5],'--','Color',darkGreen,'LineWidth',1.5);
    plot([rangeK(1)-1 rangeK(end)+1],[0.2 0.2],'--','Color',darkRed,'LineWidth',1.5);
    
    for k = 1:length(rangeK)
        for c = 1:rangeK(k)
            if effectsize(s_ind,k,c) < 0.2
                plot(rangeK(k),effectsize(s_ind,k,c),'.k','Markersize',8);
                % text(rangeK(k),effectsize(s_ind,k,c),[' ' num2str(c)]);
            end
            if effectsize(s_ind,k,c) > 0.2 && effectsize(s_ind,k,c) < 0.5
                plot(rangeK(k),effectsize(s_ind,k,c),'+','Color',darkRed,'Markersize',8);
                % text(rangeK(k),effectsize(s_ind,k,c),[' ' num2str(c)]);
            end
            if effectsize(s_ind,k,c) > 0.5 && effectsize(s_ind,k,c) < 0.8
                plot(rangeK(k),effectsize(s_ind,k,c),'o','Color',darkGreen,'Markersize',8);
                % text(rangeK(k),effectsize(s_ind,k,c),[' ' num2str(c)]);
            end
            if effectsize(s_ind,k,c) > 0.8
                plot(rangeK(k),effectsize(s_ind,k,c),'*','Color',darkBlue,'Markersize',8);
                % text(rangeK(k),effectsize(s_ind,k,c),[' ' num2str(c)]);
            end
        end
    end
    
    title([cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}],'interpreter','none')
    ylabel('Hedge''s effect size')
    xlabel('Number of FC States K')
    xticks([2 4 6 8 10 12 14 16 18 20])
    xlim([rangeK(1)-1 rangeK(end)+1])
    box off
    ylim([0 1])
    set(gca,'YTick',[0,0.2,0.5,0.8,1])
end
Save_png(Fig4, fullfile(data_dir, 'FracOccup_effetcsize.png'));
saveas(Fig4, fullfile(data_dir, 'FracOccup_effetcsize.fig'),'fig');
disp('- Plot successfully saved as FracOccup_effetcsize');

disp(' ')

end


function L = Plot_bar_pyramid(mean_P,sem_P,pval,rangeK,cond,conds)
%
% Pyramid of barplots (one row per K, one column per FC state) of the mean
% +- SEM fractional occupancy of the conditions conds, with the bars colored
% according to the p-values pval (n_K x max(K))

L = Pyramid_layout(length(rangeK),rangeK(end));
for k = 1:length(rangeK)
    for c = 1:rangeK(k)
        Pyramid_axes(L,k,c,1,[0.35 0.2]);
        hold on
        bar(1:numel(conds),mean_P(conds,k,c),'EdgeColor',Significance_color(pval(k,c),rangeK(k),rangeK),...
            'LineWidth',0.8,'FaceColor','none')
        errorbar(1:numel(conds),mean_P(conds,k,c),sem_P(conds,k,c),'LineStyle','none','Color','k',...
                 'CapSize',5,'LineWidth',0.4);
        set(gca,'XTick',1:numel(conds),'XTickLabel',[],'Fontsize',6,'TickLabelInterpreter','none','color','none')
        % Names of the conditions under the bottom row
        if k == length(rangeK)
            set(gca,'XTickLabel',cond(conds))
            xtickangle(45)
        end
        xlim([0.4 numel(conds)+0.6])
        hold off
        box off
    end
end
% Labels of the rows (K) and of the FC states (c)
Pyramid_labels(L,rangeK,1);
end
