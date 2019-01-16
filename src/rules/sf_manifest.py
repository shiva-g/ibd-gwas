"""Format sample list"""

rule format_manifest:
    input: ibd = DATA + 'raw/conrad/CAG_Data_updated12.04.2018.xlsx',
           cag1 = DATA + 'raw/cag/IRB_Microbiome_0.csv',
           cag2 = DATA + 'raw/cag/IRB_Microbiome11-2-18.csv',
           gsa_sex = DATA + 'raw/veo-ibd/veo_gender',
           gsa_hc = DATA + 'raw/veo-ibd/IJUK_GSA_719_Controls_12-1-17_PlinkFiles/IJUK_GSA_719_Controls_12-1-17.fam',
           gsa_veo = DATA + 'raw/veo-ibd/IJUK_VEOIBDx266_11-10-17PLINKFILES/IJUK_VEOIBDx266_11-10-17.fam'
    output:
        o = DATA + 'processed/MANIFEST.csv',
        d = DATA + 'processed/DISCARD_SAMPLES'
    run:
        new_cag = {'2015_CHOP_MIC_BAL_FA071_SUB':'2015_CHOP_MIC_BAL_FAM071_SUB',
                   '2015_CHOP_MIC_BAL_FAM092_sub':'2015_CHOP_MIC_BAL_FAM092_SUB',
                   '2015_CHOP_MIC_BAL_FAM098_sub':'2015_CHOP_MIC_BAL_FAM098_SUB',
                   '2015_CHOP__MIC_BAL_FAM061_SUB':'2015_CHOP_MIC_BAL_FAM061_SUB',
                   '2015_CHOP_MIC_BAL_FAM)41_SUB':'2015_CHOP_MIC_BAL_FAM041_SUB',
                   '2015_CHOP_MIC_BAL-FAM085_SUB':'2015_CHOP_MIC_BAL_FAM085_SUB',
                   '2015_CHOP_MIC_BAL_FAM097_sub':'2015_CHOP_MIC_BAL_FAM097_SUB',
                   '2915_CHOP_MIC_BAL_FAM176_SUB':'2015_CHOP_MIC_BAL_FAM176_SUB',
                   '2015_CHOP_MIC_BAL_FAM87_SUB':'2015_CHOP_MIC_BAL_FAM087_SUB',
                   '2015_CHOP_MIC_BAL_FAM94_SUB':'2015_CHOP_MIC_BAL_FAM094_SUB',
                   '2015_CHOP_MIC_BAL_FAM0104_SUB':'2015_CHOP_MIC_BAL_FAM104_SUB',
                   '2017_CHOP_BAL_VEO_007':'2018_CHOP_BAL_VEO_007',
                   '2017_CHOP_BAL_VEO_009':'2018_CHOP_BAL_VEO_009',
                   '2017_CHOP_BAL_VEO_010':'2018_CHOP_BAL_VEO_010'}
        def rename_cag(row):
            if row['SSID'] in new_cag:
                return new_cag[row['SSID']]
            return row['SSID']

        cag1 = pd.read_csv(input.cag1, skiprows=8)
        cag2 = pd.read_csv(input.cag2, skiprows=8)
        cag = pd.concat([cag1, cag2])[['Sample_ID', 'SentrixBarcode_A', 'SentrixPosition_A']].rename(columns={'Sample_ID':'SSID'})
        cag.loc[:, 'IID'] = cag.apply(lambda row: str(row['SentrixBarcode_A']) + '_' + row['SentrixPosition_A'], axis=1)
        discard_cag = ('2018_CHOP_BAL_VEO_022', '2015_CHOP_MIC_BAL_FAM106_SUB',
                       '2015_CHOP_MIC_BAL_FAM018_SUB', '2015_CHOP_MIC_BAL_FAM076_SUB',
                       '2015_CHOP_MIC_BAL_FAM101_SUB', '2015_CHOP_MIC_BAL_FAM020_SUB',
                       '2015_CHOP_MIC_BAL_FAM109_SUB'
                       )
        crit = cag.apply(lambda row: not row['SSID'] in discard_cag, axis=1)
        discard_df = cag[~crit]
        discard_df.loc[:, 'tmp'] = discard_df.apply(lambda row: '0 ' + row['IID'] + ' 0 0 0 -9', axis=1)
        discard_df[['tmp']].to_csv(output.d, index=False, header=None)
        cag.loc[:, 'SSID'] = cag.apply(rename_cag, axis=1)

        ibd = pd.read_excel(input.ibd, sheet_name='Has GWAS', skiprows=2)
        ibd.loc[:, 'race'] = ibd.apply(lambda row: 'White' if str(row['race'])=='Whtie' else str(row['race']).strip(), axis=1)
        # set sex for missing samples
        females = ['201939090006_R06C02', '201939090044_R03C01', '201939090044_R04C02', '201939090076_R03C02', '201939090076_R07C01', '201939090076_R08C02', '201939090114_R07C02',
                   '201939090156_R02C02', '201939090156_R06C02', '201939090179_R06C01', '201939090193_R04C02', '201939090193_R05C02', '201978470008_R05C02', '201939090193_R05C02']
        df = pd.merge(cag[crit], ibd, on='SSID', how='left')
        df = df.set_index('IID')
        # for sample in females:
        #     if 'nan' == str(df.at[sample, 'gender']):
        #         df.loc[sample, 'gender'] = 'Female'

        df = df.reset_index()
        df['chip'] = 'GSA+'
        def load_gsa_samples(fam_file):
            samples = {}
            with open(fam_file) as f:
                for line in f:
                    if '\t' in line:
                        sp = line.strip().split('\t')
                    else:
                        sp = line.strip().split()
                    samples[sp[1]] = True
            return samples

        gsa_hc = load_gsa_samples(input.gsa_hc)
        gsa_veo = load_gsa_samples(input.gsa_veo)
        gsa = pd.read_csv(input.gsa_sex, sep='\t', header=None, names=['IID', 'gender'])
        ls = [ [_, '-9'] for _ in gsa_hc]
        hc_df = pd.DataFrame(ls, columns=['IID', 'gender'])
        def fix_gsa_gender(row):
            assert row['gender'] in ('M', 'F'), row['gender']
            return 'Female' if row['gender']=='F' else 'Male'

        def fix_gsa_study_group(row):
            if row['IID'] in gsa_hc:
                return 'Healthy CAG'
            if row['IID'] in gsa_veo:
                return 'VEO'
            i = 1/0

        def fix_gsa_sample_type(row):
            if row['IID'] in gsa_hc:
                return 'HC'
            if row['IID'] in gsa_veo:
                return 'IBD'
            i = 1/0

        gsa.loc[:, 'gender'] = gsa.apply(fix_gsa_gender, axis=1)
        gsa = pd.concat([hc_df, gsa])
        gsa['chip'] = 'GSA'
        gsa.loc[:, 'Study Group'] = gsa.apply(fix_gsa_study_group, axis=1)
        gsa.loc[:, 'HC or IBD or ONC'] = gsa.apply(fix_gsa_sample_type, axis=1)
        pd.concat([df, gsa]).to_csv(output.o, index=False)
