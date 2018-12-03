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
        df = pd.merge(cag, ibd, on='SSID', how='left').to_csv(output.o, index=False)
