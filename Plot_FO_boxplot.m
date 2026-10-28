function Plot_FO_boxplot(stats_file,k,c,rangeK,state_color,compact)
%
% Boxplot of the fractional occupancy of one FC state in each condition, in
% the current axes, with one point per observation used in the tests
% (subject or scan), lines between the values of the same subject in
% paired comparisons of adjacent conditions, and a bar for each pair of
% conditions colored by significance (dark red p<0.05, dark green
% p<0.05/K, dark blue p<0.05/sum(rangeK)), labelled with the p-value,
% Hedges' g and its 95% bootstrap confidence interval.
%
% INPUT:
% stats_file    full path of LEiDA_Stats_FracOccup.mat
% k, c          index of K in rangeK and FC state
% rangeK        range of K
% state_color   color of the points (e.g. color of the Yeo network)
% compact       true for small panels: only the significant pairs are
%               shown, without text (default false)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

if nargin < 6 || isempty(compact)
    compact = false;
end

load(stats_file,'cond','P','P_pval2sided','Index_Conditions','Subjects','stats_level',...
     'effectsize_signed','effectsize_ci');
S = load(stats_file,'paired_test');
n_Cond = numel(cond);
[values,subjects] = FracOccup_values(P,Index_Conditions,Subjects,stats_level,n_Cond,k,c);
[cond2,cond1] = find(triu(ones(n_Cond),1)');
if isfield(S,'paired_test')
    paired = S.paired_test;
else
    paired = zeros(numel(cond1),1);
end

hold on
% Lines between the values of the same subject (paired comparisons of
% adjacent conditions)
for pair = 1:numel(cond1)
    if paired(pair) && cond2(pair) == cond1(pair) + 1
        [~,i1,i2] = intersect(subjects{cond1(pair)},subjects{cond2(pair)},'stable');
        plot([cond1(pair) cond2(pair)],[values{cond1(pair)}(i1) values{cond2(pair)}(i2)]',...
             '-','Color',[0.75 0.75 0.75],'LineWidth',0.5)
    end
end
% Boxplots and points
g = repelem(1:n_Cond,cellfun(@numel,values))';
boxplot(vertcat(values{:}),g,'Colors','k','Symbol','','Widths',0.4);
for j = 1:n_Cond
    x = j + 0.08*randn(numel(values{j}),1);
    scatter(x,values{j},18-10*compact,state_color,'filled','MarkerFaceAlpha',0.6)
end
set(gca,'XTick',1:n_Cond,'XTickLabel',cond,'FontSize',12-4*compact,'TickLabelInterpreter','none')
if ~compact
    ylabel(['Fractional occupancy (one point per ' stats_level ')'],'FontSize',12)
end

% Significance bars: p-value, Hedges' g [95% CI], colored by significance
y_max = max(vertcat(values{:}));
step = 0.12*max(y_max,eps);
n_bars = 0;
for pair = 1:numel(cond1)
    p = P_pval2sided(pair,k,c);
    if compact && p > 0.05
        continue
    end
    n_bars = n_bars + 1;
    g_val = effectsize_signed(pair,k,c);
    ci = squeeze(effectsize_ci(pair,k,c,:));
    color = Significance_color(p,rangeK(k),rangeK);
    y = y_max + n_bars*step;
    plot([cond1(pair) cond2(pair)],[y y],'-','Color',color,'LineWidth',1.2)
    if compact
        continue
    end
    if all(isfinite(ci))
        g_str = sprintf('g = %.2f [%.2f, %.2f]',g_val,ci(1),ci(2));
    else
        g_str = sprintf('g = %.2f',g_val);
    end
    if paired(pair)
        g_str = [g_str ' (paired)']; %#ok<AGROW>
    end
    text(mean([cond1(pair) cond2(pair)]),y,sprintf('p = %.2g, %s',p,g_str),'Color',color,...
         'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10)
end
ylim([0 y_max + (n_bars+1)*step])
xlim([0.4 n_Cond+0.6])
if ~compact
    title('Fractional occupancy (g > 0: left condition higher)','FontSize',11)
end
hold off
box off
