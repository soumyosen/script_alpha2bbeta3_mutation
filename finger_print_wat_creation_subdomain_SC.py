import MDAnalysis as mda
import prolif
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem

chains = ["A", "B"]
chains_A_res = ["771", "650", "596", "456"]
chains_B_res = ["63", "119", "189", "228", "329", "662"]

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

u1 = mda.Universe("last_frame_wild_noion1.pdb")
u2 = mda.Universe("last_frame_mutate_new_noion.pdb")

u1.atoms.guess_bonds()
u2.atoms.guess_bonds()

dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

dic1 = {'HBAcceptor': 'HBDonor', 'HBDonor': 'HBAcceptor', 'Cationic': 'Anionic', 'Anionic': 'Cationic', 'Hydrophobic': 'Hydrophobic', 'CationPi': 'PiCation', 'PiCation':'CationPi'}

A_dict = {range(1, 483): 'BetaProp', range(483, 633): 'Thigh', range(633, 775): 'Calf1', 
        range(775, 997): 'Calf2', range(997, 1020): 'TM', range(1020, 1039): 'CD'}

B_dict = {range(1, 83): 'PSI', range(83, 135): 'Hybrid', range(135, 379): 'BetaI', range(379, 460): 'Hybrid', 
        range(460, 462): 'PSI', range(462, 499): 'EGF1', range(499, 549): 'EGF2', range(549, 586): 'EGF3', 
        range(586, 627): 'EGF4', range(627, 632): 'Ankel', range(632, 719): 'BetaTD', range(719, 742): 'TM', 
        range(742, 788): 'CD'}


def subdomain(resnum, cha):
    new_dic = cha+'_'+'dict'
    resnum = int(resnum)
    if new_dic=="A_dict":
        subdo = { A_dict[key] for key in A_dict if resnum in key }
    else:
        subdo = { B_dict[key] for key in B_dict if resnum in key }
    subdo = " ".join(subdo)
    return subdo

def subdomain1(data_frame):
    new_dic = data_frame["chain"]+'_'+'dict'
    resnum = int(data_frame["res_num"])
    if new_dic=="A_dict":
        subdo = { A_dict[key] for key in A_dict if resnum in key }
    else:
        subdo = { B_dict[key] for key in B_dict if resnum in key }
    subdo = " ".join(subdo)
    return subdo

def column_name_change(list_of_column_names):
    new_list_of_col_AA = []
    for s in list_of_column_names:
        a1=s[0:3]
        a2=s.split('.')[0][3:]
        a3=s.split('.')[1]
        for keys, values in dic.items():
            a1 = a1.replace(keys, values)
        a4=subdomain(a2, a3)
        a1234=a3+"_"+a1+a2+"_"+a4
        new_list_of_col_AA.append(a1234)
    return new_list_of_col_AA



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
    
    ################### Backbone removal using SMART substructure
    backbone = Chem.MolFromSmarts("[C^2](=O)-[C;X4](-[H])-[N;+0]-[H]")
    backbone_pro = Chem.MolFromSmarts("[C](=O)-[C@](-[H])-[N-]")

    prot1_mol = prolif.Molecule.from_mda(prot1)
    prot1_mol = AllChem.DeleteSubstructs(prot1_mol, backbone)
    prot1_mol = AllChem.DeleteSubstructs(prot1_mol, backbone_pro)
    prot1_mol = prolif.Molecule(prot1_mol)
    lig1_mol = prolif.Molecule.from_mda(lig1)
    lig1_mol = AllChem.DeleteSubstructs(lig1_mol, backbone)
    lig1_mol = AllChem.DeleteSubstructs(lig1_mol, backbone_pro)
    lig1_mol = prolif.Molecule(lig1_mol)
    
    prot2_mol = prolif.Molecule.from_mda(prot2)
    prot2_mol = AllChem.DeleteSubstructs(prot2_mol, backbone)
    prot2_mol = AllChem.DeleteSubstructs(prot2_mol, backbone_pro)
    prot2_mol = prolif.Molecule(prot2_mol)
    lig2_mol = prolif.Molecule.from_mda(lig2)
    lig2_mol = AllChem.DeleteSubstructs(lig2_mol, backbone)
    lig2_mol = AllChem.DeleteSubstructs(lig2_mol, backbone_pro)
    lig2_mol = prolif.Molecule(lig2_mol)

    fp1 = prolif.Fingerprint(interactions="all")
    data1 = fp1.generate(lig1_mol, prot1_mol)
    data1["Frame"] = 0
    ifp1 = [data1]
    df1 = prolif.to_dataframe(ifp1, fp1.interactions.keys())

    fp2 = prolif.Fingerprint(interactions="all")
    data2 = fp2.generate(lig2_mol, prot2_mol)
    data2["Frame"] = 0
    ifp2 = [data2]
    df2 = prolif.to_dataframe(ifp2, fp2.interactions.keys())
#######################################################################################

    #fp = prolif.Fingerprint(interactions="all")
    #df1 = fp.run(u1.trajectory[:], lig1, prot1).to_dataframe()
    #df2 = fp.run(u2.trajectory[:], lig2, prot2).to_dataframe()
    #print(i)
    #print(df1.empty)
    #print(df2.empty)

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
    df["interaction1"]=df["interaction"].map(dic1)
    print(df)
    df["fp"]=df[['chain_aa', 'interaction1']].agg('_'.join, axis=1)
    df["sub_domain"] = df.apply(subdomain1, axis=1)
    df["fp_sub"]=df[["fp", "sub_domain"]].agg('_'.join, axis=1)
    #print(df)

    list_of_col_AA = []
    for s in df.columns:
        if s.isupper():
            list_of_col_AA.append(s)
    new_column_names = column_name_change(list_of_col_AA)
    

    df=df.rename(columns=dict(zip(list_of_col_AA, new_column_names)))
    nec_cols = ["fp_sub"]
    for s in new_column_names:
        nec_cols.append(s)
    #print(nec_cols) 
    df1 = df[nec_cols]
    #df1.to_csv("ft_sub_%s.csv" % i, index=False)
    #print(df1)

    df_water = pd.read_csv("water_mediate_%s.dat" % i)
    #print(df_water.empty)
    if df_water.empty == False:
        df_water=df_water.rename(columns=dict(zip(list_of_col_AA, new_column_names)))   
        df_water["sub_domain"]=df_water.apply(subdomain1, axis=1)
        df_water["amino_acid"]=df_water["amino_acid"].map(dic)
        df_water["res_num"] = df_water["res_num"].apply(str)
        df_water["aa_resnum"]=df_water["amino_acid"]+df_water["res_num"]
        df_water["chain_aa"]=df_water[['chain', 'aa_resnum']].agg('_'.join, axis=1)
        df_water["fp"]=df_water[['chain_aa', 'int_type']].agg('_'.join, axis=1)
        df_water["fp_sub"]=df_water[["fp", "sub_domain"]].agg('_'.join, axis=1)
        df_water1=df_water[nec_cols]
        df_combine = pd.concat([df1, df_water1]).reset_index(drop=True)
        df_combine.to_csv("ft_wat_sub_%s.csv" % i, index=False)
    else:
        df1.to_csv("ft_wat_sub_%s.csv" % i, index=False)












