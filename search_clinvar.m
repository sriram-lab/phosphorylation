t = readtable('variant_summary.txt');
bt = readtable('Dataset_6.csv');

unis = [];
nums = [];
%sig = [];
for i = 1:height(t)
    i
    a = split(t.Name{i});
    if length(a)>1
        if length(a{2})>7
    if strcmp(a{2}(4:6),'Tyr') & strcmp(a{2}((end-3):(end-1)),'Asp')
        %b = split(t{i,1});
        unis = [unis; t.GeneSymbol(i)];
        num = str2num(a{2}(7:(end-4)));
        if ~isempty(num)
        nums = [nums; num];
        else
            nums = [nums; 0];
        end
        %sig = [sig; a(3)];
    end
        end
    end
end

for i = 1:length(unis)
    c = find(strcmp(unis(i),bt.gene_normalized) & (nums(i)==bt.Mut_res));
    if ~isempty(c)
        bt(c,:)
        %sig(i)
    end
end

