import os
import json
try:
    import cPickle as pickle
except ImportError:
    import pickle
    
    
with open(os.getcwd() + "/uniprot_dict.json", "r") as fp:
    uniprot_dict = json.load(fp)
with open(os.getcwd() + "/uniprot2gene.json", "r") as fp:
    uniprot2gene = json.load(fp)
with open(os.getcwd() + "/corona_id_to_seq.json", "r") as fp:
    corona_id_to_seq = json.load(fp)
with open(os.getcwd() + "/corona_id_to_protein.json", "r") as fp:
    corona_id_to_protein = json.load(fp)
with open(os.getcwd() + "/organism_to_keys.json", "r") as fp:
    organism_to_keys = json.load(fp)
    
    
def longest_substring(s1, s2):
    t = [[0]*(1+len(s2)) for i in range(1+len(s1))]
    l, xl = 0, 0
    for x in range(1,1+len(s1)):
        for y in range(1,1+len(s2)):
            if s1[x-1] == s2[y-1]:
                t[x][y] = t[x-1][y-1] + 1
                if t[x][y]>l:
                    l = t[x][y]
                    xl  = x
            else:
                t[x][y] = 0
    return s1[xl-l: xl]


def find_viral_species(key_to_find):
    for key, value in organism_to_keys.items():
        if key_to_find in value:
            return key

def build_dump():
    temp_l = []
    for corona_id, corona_seq in corona_id_to_seq.items():
        print(corona_id)
        for human_id, human_seq in uniprot_dict.items():
            match = longest_substring(corona_seq, human_seq)
            if match:
                temp_dict = dict()
                temp_dict["peptide"] = match
                temp_dict["peptide_length"] = len(match)
                temp_dict["human_protein"] = uniprot2gene[human_id]
                temp_dict["human_protein_annotation"] = human_id
                temp_dict["viral_species"] = find_viral_species(corona_id)
                temp_dict["viral_protein"] = corona_id_to_protein[corona_id]
                temp_dict["viral_protein_annotation"] = corona_id
                temp_dict["human_protein_peptide_start"] = uniprot_dict[human_id].find(match)
                temp_dict["human_protein_peptide_end"] = temp_dict["human_protein_peptide_start"] + temp_dict["peptide_length"] -1
                temp_dict["viral_protein_peptide_start"] = corona_id_to_seq[corona_id].find(match)
                temp_dict["viral_protein_peptide_end"] = temp_dict["viral_protein_peptide_start"] + temp_dict["peptide_length"] -1
                temp_l.append(temp_dict)
    return temp_l

dump = build_dump()

with open(os.getcwd() + "/corona_dump.p", "wb") as fp:
    pickle.dump(temp_l, fp)