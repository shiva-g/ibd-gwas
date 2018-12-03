rule format_manifest:
    input: ibd = DATA + 'raw/conrad/CAGData112818.xlsx',
           cag1 = DATA + 'raw/cag/IRB_Microbiome_0.csv',
           cag2 = DATA + 'raw/cag/IRB_Microbiome11-2-18.csv'
    output:
        o = DATA + 'processed/MANIFEST.csv'
    run:
        ibd = pd.read_excel(input.ibd, sheet_name='Has GWAS', skiprows=2)
        cag1 = pd.read_csv(input.cag1, skiprows=8)
        cag2 = pd.read_csv(input.cag2, skiprows=8)
        cag = pd.concat([cag1, cag2])[['Sample_ID', 'SentrixBarcode_A', 'SentrixPosition_A']].rename(columns={'Sample_ID':'SSID'})
        cag.loc[:, 'IID'] = cag.apply(lambda row: str(row['SentrixBarcode_A']) + '_' + row['SentrixPosition_A'], axis=1)
        ibd.loc[:, 'race'] = ibd.apply(lambda row: 'White' if str(row['race'])=='Whtie' else str(row['race']).strip(), axis=1)
        # set sex for missing samples
        females = ['201939090006_R06C02', '201939090044_R03C01', '201939090044_R04C02', '201939090076_R03C02', '201939090076_R07C01', '201939090076_R08C02', '201939090114_R07C02',
                   '201939090156_R02C02', '201939090156_R06C02', '201939090179_R06C01', '201939090193_R04C02', '201939090193_R05C02', '201978470008_R05C02', '201939090193_R05C02']
        df = pd.merge(cag, ibd, on='SSID', how='left')
        df = df.set_index('IID')
        for sample in females:
            if 'nan' == str(df.at[sample, 'gender']):
                df.loc[sample, 'gender'] = 'Female'

        df.reset_index().to_csv(output.o, index=False)
