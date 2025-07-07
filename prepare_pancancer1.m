t = readtable('mmc1.xlsx','Sheet','Table S1C');

idxes = [];
for i = 1:height(t)
    i/height(t)
    if ~isempty(t{i,20}{:})
        a = split(t{i,20},{'[',']'});
        a = split(a{2},'|');
        idx1 = [];
        for j = 1:length(a)
            idx = find(isstrprop(a{j},"alpha"),1)-1*(j==1);
            idx1 = [idx1, num2str(idx),','];
        end
        idxes = [idxes; {idx1}];
    else
        idxes = [idxes; {0}];
    end
end

t2 = readtable('mmc2.xlsx','Sheet','Table S2D');
phosrows = strcmp(t2{:,6},'phosphoproteome') | strcmp(t2{:,6},'ptm');
np_to_get = t2(phosrows,1);
uni_to_get = t2(phosrows,3);
np0 = t{:,1};
np = {};
for i = 1:length(np0)
    i/length(np0)
    np1 = split(np0{i},'_');
    np1 = [np1{1}, '_', np1{2}];
    np = [np; np1];
end

for_table = [];
for i = 1:height(np_to_get)
    i/height(np_to_get)
    a = t{find(strcmp(np_to_get{i,:},np)),[18,20]};
    b = idxes(find(strcmp(np_to_get{i,:},np)));
    for_table = [for_table; [repmat(np_to_get{i,:},size(a,1),1), repmat(uni_to_get{i,:},size(a,1),1), a, b]];
end
    
writecell(for_table,'phostable.csv')
writecell(unique(for_table(:,2)),'phostable_unis.csv')