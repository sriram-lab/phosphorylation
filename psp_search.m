t = readtable('~/Desktop/cancer_ptm/alpha_predictions_1.csv'); %set to path of alpha_predictions
tp = readtable('phosite');
resis = split(tp{:,5},'-');
resis = resis(:,1);
restype = cell(1,length(resis));
resnum = zeros(1,length(resis));
for i = 1:length(resis)
    restype(i) = {resis{i}(1)};
    resnum(i) = str2num(resis{i}(2:end));
end
organism = tp{:,7};
id = tp{:,3};

id_sel = id(strcmp(organism,'human') & strcmp(restype,'Y')');
num_sel = resnum(strcmp(organism,'human') & strcmp(restype,'Y')');

all_id = split(t.pdb_name,'-');
all_id = all_id(:,2);
all_num = t.resinum;

xgb = [];
unis = [];
for i = 1:length(id_sel)
    i/length(id_sel)
    j = find(strcmp(id_sel(i),all_id) & num_sel(i)==all_num);
    if ~isempty(t.XGB_Predictions(j))
    unis = [unis;all_id(i)];
    xgb = [xgb; t.XGB_Predictions(j(1))];
    end
end
