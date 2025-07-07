
phoslist = [];
for_table = readcell('phostable.csv');
match = -1*ones(1,length(for_table));

for i = 1:length(for_table)
    i/length(for_table)
    
    uni = for_table{i,2};
    aa = for_table{i,3};
    if length(aa) > 0
        aa = split(aa,{'[',',',']'});
        aa = aa{2};
        aa = split(aa,' ');
        
        for j = 1:length(aa)
            aaj = aa{j};
            d = isstrprop(aaj,"alphanum");
            aaj = aaj(d);
            if isfile(sprintf('./seqalign/%s_%s.csv',uni,aaj))
                seq2 = readcell(sprintf('./seqalign/%s_%s.csv',uni,aaj));
                seq2 = seq2{2};
                
                if for_table{i,5}~=0
                idxes = split(for_table{i,5},',');
                if j <= length(idxes)
                if ~isempty(idxes{j})
                    b = find(isstrprop(seq2,"alpha"),9-str2num(idxes{j}));
                    match(i) = str2num(aaj(isstrprop(aaj,"digit")))==b(end);
                    if aaj(1)=='S' | aaj(1)=='T' | aaj(1)=='Y'
                    phoslist = [phoslist; {num2str(str2num(aaj(isstrprop(aaj,"digit")))==b(end))},uni,aaj(1),num2str(b(end)),aaj(2:(end-1))];
                    end
                end
                end
                end
            end
        end
    end
end

a = table(phoslist(:,1),phoslist(:,2),phoslist(:,3),phoslist(:,4),phoslist(:,5));

writetable(unique(a),'phoslist_with_conversion.csv')