"""Derive IBD relevant pathways
"""
import pandas as pd
import argparse

def is_redox(row):
    if 'NF_KAPPAB_SIGNALING' in row['GO'] or 'OXIDATIVE_STRESS' in row['GO']:
        return True
    return False

def is_er_stress(row):
    if 'ENDOPLASMIC_RETICULUM_STRESS' in row['GO']:
        return True
    for gene in ('CPEB4', 'ORMDL3', 'SERINC3', 'XBP1', ):
        if gene in row['gene_ls']:
            return True
    return False

def is_autophagy(row):
    if 'GO_AUTOPHAGY' in row['GO']:
        return True
    for gene in ('IRGM', 'NOD2', 'ATG16L1', 'LRPK2', 'CUL2', 'PARK7', 'DAP'):
        if gene in row['gene_ls']:
            return True
    return False

def is_apoptosis(row):
    if 'APOPTOTIC_SIGNALING_PATHWAY' in row['GO']:
        return True
    for gene in ('FASLG', 'THADA', 'DAP', 'PUS10', 'MST1', ):
        if gene in row['gene_ls']:
            return True
    return False

def is_innate_immune(row):
    if 'GO_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE' in row['GO']:
        return True
    return False

def is_cell_migration(row):
    if 'CELL_MIGRATION' in row['GO']:
        return True
    for gene in ('ARPC2', 'LSP1', 'AAMP', ):
        if gene in row['gene_ls']:
            return True
    return False

def is_tcell_regulation(row):
    if 'T_CELL' in row['GO']:
        return True
    for gene in ('NDFIP1', 'TNFSF8', 'TAGAP', 'IL2', 'IL2R', 'TNRFSF9', 'PIM3', 'IL7R', 'IL12B', 'PRDM1', 'ICOSLG', 'TNFSF8', 'IFNG'):
        if gene in row['gene_ls']:
            return True
    return False

def is_bcell_regulation(row):
    if 'B_CELL' in row['GO']:
        return True
    for gene in ('IL5', 'IKZF1', 'BACH2', 'IL7R', 'IRF5', ):
        if gene in row['gene_ls']:
            return True
    return False

def mk_pathways(row):
    pathways = {'er_stress':is_er_stress,
                'bcell_regulation':is_bcell_regulation,
                'tcell_regulation':is_tcell_regulation,
                'apoptosis':is_apoptosis,
                'cell_migration':is_cell_migration,
                'redox':is_redox,
                'autophagy':is_autophagy,
                'innate_immunity':is_innate_immune}
    paths = ','.join([p for p in pathways if pathways[p](row)])
    if not paths:
        paths = 'none'
    return paths 

def mk_genes(row):
    genes = []
    if 'ANN=' in str(row['eff']):
        ls = row['eff'].split('ANN=')[1].split(';')[0].split(',')
        genes = set([x.split('|')[3] for x in ls])
    return genes

def main(args):
    df = pd.read_csv(args.assoc_table, sep='\t')
    df.loc[:, 'gene_ls'] = df.apply(mk_genes, axis=1)
    df.loc[:, 'pathways'] = df.apply(mk_pathways, axis=1)
    df.to_csv(args.out, index=False, sep='\t')

if __name__ == "__main__":
    desc = 'Derive pathways from GO and papers.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('assoc_table', 'out',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
