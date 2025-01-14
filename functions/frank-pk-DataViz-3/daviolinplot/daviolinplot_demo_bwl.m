% daviolinplot_demo a few examples of daviolinplot functionality 
%
% Povilas Karvelis
% 05/05/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
rng('default')

% data in a matrix (+ grouping indices)
data2 = [randn([30,4]); randn([30,4]);...
         randn([30,4]); randn([30,4])];
group_inx = [ones(1,30), 2.*ones(1,30) 3.*ones(1,30) 4.*ones(1,30)];

% group_names = {'Humans', 'Dogs' , 'God', 'Potato'};
% condition_names = {'Water', 'Land', 'Moon', 'Hyperspace'};


load('within_n_between_subj_correlation.mat');

data2 = [betweenCorr{1,1}.rest, betweenCorr{1,1}.zerobk, betweenCorr{1,1}.twobk];
data2 = [data2; withinCorr{1,1}.rest, withinCorr{1,1}.zerobk, withinCorr{1,1}.twobk];

group_inx=[ones(1,length(betweenCorr{1,1}.rest)), 2.*ones(1,length(withinCorr{1,1}.rest))];



group_names = {'Between', 'Within'};
condition_names = {'Rest', 'Zerobk', 'Twobk'};


% an alternative color scheme for some plots
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 

figure;
h = daviolinplot(data2(:,1:3),'groups',group_inx, 'outliers', 0, 'outsymbol','k+',...
    'xtlabels', condition_names,'color',c,'scatter',2,'jitter',1,...
    'box',1,'boxcolors','same','scattercolors','same',...
    'boxspacing',1.1,'legend',group_names(1:2));
ylabel('Performance');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',10);





