t = readtable('Dataset_6.csv');
%a = t.Area;
%b = t.phosphorylation;
%corr(a(~isnan(a) & ~isnan(b)),b(~isnan(a) & ~isnan(b)))

unis = [];
genes = [];
t3 = [];
t2 = readtable('mmc2_2.xlsx','sheet','Table S1S');
for i = 1:height(t)
    n = find(strcmp(t.gene_normalized{i},t2{:,1}));
    unis = [unis; t.uniprot(i)];
    t3 = [t3; t.predictions(i), mean(t2.degree(n)), mean(t2.betweenness(n)), mean(t2.closeness(n)), mean(t2.pagerank(n)), mean(t2.geneKO(n))];
end
unis(isnan(t3(:,2)),:) = [];
t3(isnan(t3(:,2)),:) = [];
close all
figure(1)
boxplot(t3(:,2),t3(:,1)>1)
title('degree')
figure(2)
boxplot(log(t3(:,3)),t3(:,1)>1)
title('betweenness')
figure(3)
boxplot(t3(:,4),t3(:,1)>1)
title('closeness')
figure(4)
boxplot(t3(:,5),t3(:,1)>1)
title('pagerank')
figure(5)
boxplot(t3(:,6),t3(:,1)>1)
title('Gene knockout')

for i = 1:5
    figure(i)
    ylabel('Centrality')
    xlabel('\Delta\DeltaG > 1 kcal/mol')
    set(gca,'xticklabel',{'no','yes'})
    set(gca,'Fontsize',20)
end

[h,p_degree] = ttest2(t3(t3(:,1)>1,2),t3(t3(:,1)<=1,2)) % degree
a = log(t3(:,3));
[h,p_betweenness] = ttest2(a(~isnan(a) & (t3(:,1)>1) & ~isinf(a)),a(~isnan(a) & (t3(:,1)<=1)& ~isinf(a))) % betweenness
[h,p_closeness] = ttest2(t3(t3(:,1)>1,4),t3(t3(:,1)<=1,4)) % closeness
[h,p_pagerank] = ttest2(t3(t3(:,1)>1,5),t3(t3(:,1)<=1,5)) % pagerank
[h,p_KO] = ttest2(t3(t3(:,1)>1,6),t3(t3(:,1)<=1,6)) % geneKO
%histogram(t3(t3(:,1)>1,2),0:0.1:1,'Normalization','Probability')
%hold on
%histogram(t3(t3(:,1)<=1,2),0:0.1:1,'Normalization','Probability')
u = unis((t3(:,1)>1) & (t3(:,2)<200));
genes = [];
for i = 1:length(u)
    genes = [genes; t.gene_normalized(strcmp(u{i},t{:,2}))];
end

t_essential = readtable('CSEGs_CEGs.txt');
essential_genes = t_essential.gene(strcmp(t_essential.essentiality,'CEGs'));
essential_genes_cancer = t_essential.gene(strcmp(t_essential.essentiality,'CSEGs'));

eg_column = zeros(1,height(t));
for i = 1:height(t)
    if ~isempty(t.gene_normalized(i))
        if sum(strcmp(t.gene_normalized(i),essential_genes))
            eg_column(i) = 1;
        end
    end
end
  
boxplot(t.predictions,eg_column)
[h,p] = ttest2(t.predictions(eg_column==0),t.predictions(eg_column==1))

for i = 1:5
figure(i)
set(gca,'box','off')
set(gca,'linewidth',1.3)
set(gcf,'color','white')
grid('on')
    set(gca,'fontsize',18)
set(findobj(gca,'Type','text'),'FontSize',18)
set(findobj(gca,'Type','line'),'linewidth',1.75,'color','k')
 set(gca,'TickDir', 'out')
set(gca,'fontname','helvetica','fontweight','normal')
end
      
%exportgraphics(gcf, 'bp_degree.tif','Resolution',300)