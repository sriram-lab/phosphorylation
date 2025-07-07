t = readtable('Dataset_6.csv');
todel = zeros(1,height(t));
for i = 2:height(t)
    if strcmp(t.uniprot(i),t.uniprot(i-1)) && t.Mut_res(i)==t.Mut_res(i-1)
        todel(i) = 1;
    end
end
todel(isnan(tpredictions))=1;
t(todel==1,:) = [];
writetable(t,'Dataset_6.csv');