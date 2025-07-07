tsu = readtable('dG_non_redundant_natural_Fig5.csv');

ymut = strcmp(tsu{:,3},'Y');
ydis = -tsu{ymut,5} + tsu{ymut,20};
%ypos = tsu{ymut,2};
%yaa = tsu{ymut,3};
%yfile = tsu{ymut,1};

ymut2 = strcmp(tsu{:,3},'Y');
ydis2 = -tsu{ymut2,21} + tsu{ymut2,20};
%ypos2 = tsu{ymut2,2};
%yaa2 = tsu{ymut2,3};
%yfile2 = tsu{ymut2,1};

smut = strcmp(tsu{:,3},'S');
sdis = -tsu{smut,5} + tsu{smut,12};
%spos = tsu{smut,2};
%saa = tsu{smut,3};
%sfile = tsu{smut,1};

tmut = strcmp(tsu{:,3},'T');
tdis = -tsu{tmut,5} + tsu{tmut,11};
%tpos = tsu{tmut,2};
%taa = tsu{tmut,3};
%tfile = tsu{tmut,1};

%histogram(sdis,-3:0.5:5,'Normalization','Probability')
%hold on
%histogram(tdis,-3:0.5:5,'Normalization','Probability')
%histogram(ydis,-3:0.5:5,'Normalization','Probability')

boxplot([sdis; tdis; ydis; ydis2],[ones(length(sdis),1);2*ones(length(tdis),1); 3*ones(length(ydis),1);4*ones(length(ydis),1)]);

writetable(table([tsu{smut,3}; tsu{tmut,3}; tsu{ymut,3}],[sdis; tdis; ydis]),'tsutable.csv')
writetable(table([sdis; tdis; ydis],[sfile; tfile; yfile],[saa; taa; yaa],[spos; tpos; ypos]),'phoslist_tsu.csv')
%updb = unique([sfile; tfile; yfile]);
%for i = 1:length(updb)
%    updb{i} = updb{i}(1:4);
%end
%updb
%writecell(updb,'tsu_pdbs.csv')