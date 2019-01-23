"""Format sample list"""

rule format_manifest_gsa:
    input:
        gsa_hc = DATA + 'raw/veo-ibd/IJUK_GSA_719_Controls_12-1-17_PlinkFiles/IJUK_GSA_719_Controls_12-1-17.fam',
        gsa_veo = DATA + 'raw/veo-ibd/IJUK_VEOIBDx266_11-10-17PLINKFILES/IJUK_VEOIBDx266_11-10-17.fam',
        veo_dat = DATA + 'raw/veo-ibd/SNPList_01162019.xls'
    output:
        o = DATA + 'interim/manifest/gsa.csv',
    run:
        def fix_race(row):
            if row['race'] in ('White', 'white'):
                return 'White'
            if row['race'] in ('Black / African American', 'Black/African American'):
                return 'Black/African American'
            if row['race'] in ('UN', 'Declined to Answer', 'Declined to answer') or str(row['race']).strip() == '' or str(row['race'])=='nan':
                return 'Unknown'
            if row['race'] in ('White,Black/African American', 'White,Black / African American', 'Black / African American,White'):
                return 'White, Black/African American'
            return row['race']

        def fix_veo_id(row):
            """Add 0s"""
            assert 'A' in str(row['IID'])
            id = row['IID'].strip('A')
            l = len(id)
            return '0'*l + id

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

        def check_veo(samples, veo_df, dup_df):
            """samples are from the veo fam file
              veo_df and dup_df are from the
              sample data xls
            """
            xls_samples = set(veo_df['IID'].values) | set(dup_df['IID'].values)
            fam_samples = set(samples)
            not_in_xls = fam_samples - xls_samples
            # add assert back when learn about 6152A
            #assert not len(not_in_xls), not_in_xls
            not_in_fam = xls_samples - fam_samples
            assert not len(not_in_fam), not_in_fam

        def update_snp_xls_id(row):
            if row['IID'] in ('0739A', '0295A', '0800A'):
                return row['IID']
            else:
                return row['IID'].lstrip('0')

        veo_df = pd.read_excel(input.veo_dat, sheet_name='SNP Array').rename(columns={'Subject ID':'IID', 'Race':'race', 'Ethnicity':'ethnicity', 'Sex':'gender'})[['IID', 'race', 'ethnicity', 'gender']]
        veo_df.loc[:, 'IID'] = veo_df.apply(update_snp_xls_id, axis=1)
        veo_df['Study Group'] = 'VEO'
        veo_df['HC or IBD or ONC'] = 'IBD'
        dup_df = pd.read_excel(input.veo_dat, sheet_name='Patient duplicates').rename(columns={'Subject ID':'IID'})
        dup_df.loc[:, 'IID'] = dup_df.apply(update_snp_xls_id, axis=1)
        gsa_dups = dup_df['IID'].values

        gsa_hc = load_gsa_samples(input.gsa_hc)
        gsa_veo = load_gsa_samples(input.gsa_veo)
        check_veo(gsa_veo, veo_df, dup_df)

        ls = [ [_, 'Unknown', 'Unkown', '-9', 'Healthy CAG', 'HC'] for _ in gsa_hc]
        hc_df = pd.DataFrame(ls, columns=['IID', 'race', 'ethnicity', 'gender', 'Study Group', 'HC or IBD or ONC'])

        df = pd.concat([hc_df, veo_df])
        df['chip'] = 'GSA'
        df.loc[:, 'is_gsa_duplicate'] = df.apply(lambda row: row['IID'] in gsa_dups, axis=1)
        df.loc[:, 'race'] = df.apply(fix_race, axis=1)
        df.to_csv(output.o, index=False, sep='\t')

rule format_manifest_gsaplus:
    input: ibd = DATA + 'raw/conrad/CAG_Data_updated12.04.2018.xlsx',
           cag1 = DATA + 'raw/cag/IRB_Microbiome_0.csv',
           cag2 = DATA + 'raw/cag/IRB_Microbiome11-2-18.csv',
    output:
        o = DATA + 'interim/manifest/gsaplus.csv',
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

        def fix_race(row):
            if row['race'] in ('White', 'white'):
                return 'White'
            if row['race'] in ('Black / African American', 'Black/African American'):
                return 'Black/African American'
            if row['race'] in ('UN', 'Declined to Answer', 'Declined to answer') or str(row['race']).strip() == '' or str(row['race'])=='nan':
                return 'Unknown'
            if row['race'] in ('White,Black/African American', 'White,Black / African American', 'Black / African American,White'):
                return 'White, Black/African American'
            return row['race']

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

        ibd = pd.read_excel(input.ibd, sheet_name='Has GWAS', skiprows=2).rename(columns={'study_id':'SubjectID'})
        ibd.loc[:, 'race'] = ibd.apply(lambda row: 'White' if str(row['race'])=='Whtie' else str(row['race']).strip(), axis=1)
        # set sex for missing samples
        females = ['201939090006_R06C02', '201939090044_R03C01', '201939090044_R04C02', '201939090076_R03C02', '201939090076_R07C01', '201939090076_R08C02', '201939090114_R07C02',
                   '201939090156_R02C02', '201939090156_R06C02', '201939090179_R06C01', '201939090193_R04C02', '201939090193_R05C02', '201978470008_R05C02', '201939090193_R05C02']
        df = pd.merge(cag[crit], ibd, on='SSID', how='left')
        df = df.set_index('IID')

        df = df.reset_index()
        df['chip'] = 'GSA+'

        df.loc[:, 'is_gsa_duplicate'] = False
        df.loc[:, 'race'] = df.apply(fix_race, axis=1)
        df.to_csv(output.o, sep='\t', index=False)

rule format_manifest:
    input:
        expand(DATA + 'interim/manifest/{chip}.csv', chip=('gsa', 'gsaplus'))
    output:
        o = DATA + 'processed/MANIFEST.csv'
    run:
        pd.concat([pd.read_csv(i, sep='\t') for i in input]).to_csv(output.o, index=False)
