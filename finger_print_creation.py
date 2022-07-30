import MDAnalysis as mda
import prolif

chains = ["A", "B"]
chains_A_res = ["1032", "771", "650", "596", "63"]
chains_B_res = ["189", "433", "654", "662"]

strings = []

for i in chains_A_res:
    string = "segid %s and resid %s" % (chains[0], i)
    strings.append(string)

for i in chains_B_res:
    string = "segid %s and resid %s" % (chains[1], i)
    strings.append(string)

#print(strings)

#new_string = "same residue as protein and not (%s) and around 10.0 (%s)" % (strings[0], strings[0])
#print(new_string)

u1 = mda.Universe("a2bb3_str_Hadded_min.pdb")
u2 = mda.Universe("a2bb3_str_mutated_Hadded_min.pdb")

u1.atoms.guess_bonds()
u2.atoms.guess_bonds()

dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


for i in range(len(strings)):
#for i in range(1):
    #print(strings[i])
    prot1 = u1.select_atoms("same residue as protein and not (%s) and around 10.0 (%s)" % (strings[i], strings[i]))
    lig1 = u1.select_atoms("%s" % strings[i])
     
    prot2 = u2.select_atoms("same residue as protein and not (%s) and around 10.0 (%s)" % (strings[i], strings[i]))
    lig2 = u2.select_atoms("%s" % strings[i])

    #print(prot1)
    #print(lig1)
    #print(prot2)
    #print(lig2)

    fp = prolif.Fingerprint(interactions="all")
    df1 = fp.run(u1.trajectory[:], lig1, prot1).to_dataframe()
    df2 = fp.run(u2.trajectory[:], lig2, prot2).to_dataframe()

    df1 = df1.unstack(level=-1).reset_index()
    df1 = df1.drop("Frame", axis=1)
    df1 = df1.rename(columns={0: "Present"})
    df1 = df1.pivot(index=["protein", "interaction"], columns="ligand", values="Present").reset_index()

    df2 = df2.unstack(level=-1).reset_index()
    df2 = df2.drop("Frame", axis=1)
    df2 = df2.rename(columns={0: "Present"})
    df2 = df2.pivot(index=["protein", "interaction"], columns="ligand", values="Present").reset_index()

    df = df1.merge(df2, how="outer", on=["protein", "interaction"])
    df = df.fillna(False)
    df["amino_acid"]=df["protein"].str[0:3]
    df["res_num"]=df["protein"].str.split('.').str[0].str[3:]
    df["chain"]=df["protein"].str.split('.').str[1]
    df["amino_acid"]=df["amino_acid"].map(dic)
    df["amino_acid"]=df["amino_acid"]+df["res_num"]
    df["chain_aa"]=df[['chain', 'amino_acid']].agg('_'.join, axis=1)
    df["fp"]=df[['chain_aa', 'interaction']].agg('_'.join, axis=1)
    #print(df)

    list_of_col_AA = []
    for s in df.columns:
        if s.isupper():
            list_of_col_AA.append(s)
    #print(list_of_col_AA)

    new_list_of_col_AA = []
    for s in list_of_col_AA:
        a1=s[0:3]
        a2=s.split('.')[0][3:]
        a3=s.split('.')[1]
        for keys, values in dic.items():
            a1 = a1.replace(keys, values)
        a123=a3+"_"+a1+a2
        new_list_of_col_AA.append(a123)
    #print(new_list_of_col_AA)

    df=df.rename(columns=dict(zip(list_of_col_AA, new_list_of_col_AA)))
    nec_cols = ["fp"]
    for s in new_list_of_col_AA:
        nec_cols.append(s)
    #print(nec_cols) 
    df1 = df[nec_cols]
    df1.to_csv("ft_%s.csv" % i, index=False)
    #print(df1)















