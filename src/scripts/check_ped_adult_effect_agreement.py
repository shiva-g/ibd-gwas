"""Flip ped alleles and OR to match adult.
   Flag effect agreement.
"""
import pandas as pd
import argparse

def mk_ped_alleles(row, ped_a1, ped_a2, ped_or):
    a1, a2, oR = '', '', ''
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    nuc1, nuc2 = row[ped_a1], row[ped_a2]
    if row['A1_adult'] == nuc1 and row['A2_adult'] == nuc2:
        a1, a2, oR = nuc1, nuc2, row[ped_or]
    elif row['A1_adult'] == nuc2 and row['A2_adult'] == nuc1:
        a1, a2, oR = nuc2, nuc1, 1/row[ped_or]

    if not a1:
        # test complement
        nuc1, nuc2 = comp[row[ped_a1]], comp[row[ped_a2]]
        if row['A1_adult'] == nuc1 and row['A2_adult'] == nuc2:
            a1, a2, oR = nuc1, nuc2, row[ped_or]
        elif row['A1_adult'] == nuc2 and row['A2_adult'] == nuc1:
            a1, a2, oR = nuc2, nuc1, 1/row[ped_or]
        else:
            print(row)
            print(row['A1_adult'], nuc1, row['A2_adult'], nuc2)
            i =1/0

    dat = {'a1':a1, 'a2':a2, 'or':oR}
    return pd.Series(dat, index=['a1', 'a2', 'or'])

def or_agree(row):
    """Does OR agree between pediatric snptest/plink and adult?"""
    if row['OR_snptest'] > 1 and row['OR_plink'] < 1:
        return 'strange plink/snptest mismatch'
        #i = 1/0
    if row['OR_snptest'] < 1 and row['OR_plink'] > 1:
        return 'strange plink/snptest mismatch'
        #print(row)
        #i = 1/0

    if row['OR_snptest'] > 1 and row['OR_adult'] > 1:
        return 'match'
    if row['OR_snptest'] < 1 and row['OR_adult'] < 1:
        return 'match'

    return 'mismatch'

def main(args):
    df = pd.read_csv(args.assoc_table, sep='\t').dropna()
    df.loc[:, 'OR_snptest'] = df.apply(lambda row: 1/row['OR_snptest'], axis=1)

    df[['a1', 'a2', 'or']] = df.apply(lambda row: mk_ped_alleles(row, 'A1_plink', 'A2_plink', 'OR_plink'), axis=1)
    df.loc[:, 'A1_plink'] = df['a1']
    df.loc[:, 'A2_plink'] = df['a2']
    df.loc[:, 'OR_plink'] = df['or']
    df = df.drop(['a1', 'a2', 'or'], axis=1)

    df[['a1', 'a2', 'or']] = df.apply(lambda row: mk_ped_alleles(row, 'A1_snptest', 'A2_snptest', 'OR_snptest'), axis=1)
    df.loc[:, 'A1_snptest'] = df['a1']
    df.loc[:, 'A2_snptest'] = df['a2']
    df.loc[:, 'OR_snptest'] = df['or']
    df = df.drop(['a1', 'a2', 'or'], axis=1)
    df.loc[:, 'OR_agree'] = df.apply(or_agree, axis=1)
    df.to_csv(args.out, index=False, sep=',')

if __name__ == "__main__":
    desc = 'Flip ped alleles and OR to agree w/ adult, and flag ped/adult agreement.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('assoc_table', 'out',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
