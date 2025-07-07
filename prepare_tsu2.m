tsu = readtable('../dG_non_redundant_natural_Fig5.csv');
updbs = unique(tsu{:,1});

tsu_seqs = [];
for i = 1:length(updbs)
    sel = strcmp(updbs{i},tsu{:,1});
    aas_cell = tsu{sel,3}';
    aas = '';
    for i = 1:length(aas_cell)
        aas = [aas, aas_cell{i}];
    end
    tsu_seqs = [tsu_seqs;{aas}];
end

map_from = {'ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','SEC','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP'};
map_to = {'R','H','K','D','E','S','T','N','Q','C','U','G','P','A','V','I','L','M','F','Y','W'};
m = containers.Map(map_from,map_to);

pdb_seqs = [];
pdb_resis = [];
for i = 1:length(updbs)
    i
    seq = '';
    prot = pdbread(sprintf('../pdbfiles2/%s',updbs{i}));
    resis = [prot.Model.Atom.resSeq];
    uresis = unique(resis);
    pdb_resis = [pdb_resis; {uresis}];
    for j = 1:length(uresis)
        n = find(uresis(j)==resis,1);
        seq = [seq, m(prot.Model.Atom(n).resName)];
    end
    pdb_seqs = [pdb_seqs; {seq}];
end

aln_struct = struct;
for i = 1:length(pdb_seqs)
    [a,aln] = nwalign(tsu_seqs{i},pdb_seqs{i});
    aln_struct(i).alignment = aln;
end


t = readtable('./phoslist_tsu.csv');

resnum_pdb = [];
for j = 1:height(t)
    pdbnum = find(strcmp(t{j,2},updbs));
    resnum_tsu = t{j,4};
    
    i = 1;
    cnt = 0;
    while cnt < resnum_tsu
        if ~strcmp(aln_struct(pdbnum).alignment(1,i),'-')
            cnt = cnt + 1;
        end
        i = i + 1;
    end
    
    if i - sum(aln_struct(pdbnum).alignment(3,1:(i-1))=='-')-1 > 0
        if (strcmp(t{j,3},pdb_seqs{pdbnum}(i - sum(aln_struct(pdbnum).alignment(3,1:(i-1))=='-')-1)))
            resnum_pdb = [resnum_pdb; pdb_resis{pdbnum}(i - sum(aln_struct(pdbnum).alignment(3,1:(i-1))=='-')-1)];
        else
            pdbnum, i
            aln_struct(pdbnum).alignment(:,i-1)
            resnum_pdb = [resnum_pdb; 0];
        end
    else
        resnum_pdb = [resnum_pdb; 0];
    end
end



t_new = [t, table(resnum_pdb)];
t_new(find(resnum_pdb==0),:) = [];
writetable(t_new,'phoslist_tsu_new.csv')