t = readtable('~/Desktop/cancer_ptm/alpha_predictions_1.csv'); %set to path of alpha_predictions
t2 = readtable('Dataset_6.csv');
unis = split(t.pdb_name,'-');
unis = unis(:,2);
newcol = zeros(1,height(t2));
for i = 1:height(t2)
    uni = t2.uniprot(i);
    res = t2.Mut_res(i);
    index = find(strcmp(uni,unis) & (t.resinum==res));
    if isempty(t.XGB_Predictions(index))
        i
    else
        newcol(i) = t.XGB_Predictions(index);
    end
end

%newcol = table(newcol','VariableNames',{'rapid_method'});
%newtable = [t2,newcol];
%writetable(newtable,'final_output_incl_rapid.csv')
