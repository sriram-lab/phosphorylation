t = readtable('humsavar2.txt','ReadVariableNames',false);
bt = readtable('Dataset_6.csv');

unis = [];
nums = [];
sig = [];
for i = 1:height(t)
    a = split(t.Var2{i});
    if length(a)>1
    if strcmp(a{2}(3:5),'Tyr') & strcmp(a{2}((end-2):end),'Asp')
        b = split(t{i,1});
        unis = [unis; b(2)];
        nums = [nums; str2num(a{2}(6:(end-3)))];
        sig = [sig; a(3)];
    end
    end
end

for i = 1:length(unis)
    c = find(strcmp(unis(i),bt.UniprotID) & (nums(i)==bt.ResNum));
    if ~isempty(c)
        bt(c,:)
        sig(i)
    end
end

