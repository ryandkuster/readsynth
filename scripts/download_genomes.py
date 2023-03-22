import os
import pandas as pd
import subprocess
import sys


'''
this script requires an up-to-date 'assembly_summary_genbank.txt'
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

sys.argv[1] is assembly_summary_genbank.txt file
sys.argv[2] is a tab_separated file from bracken output
sys.argv[3] is the location to store the downloaded genome
'''


def main():
    with open('summary.txt', 'w') as f1, open('abundances.csv', 'w') as f2:
        if not os.path.exists(sys.argv[3]):
            sys.exit('outdir does not exist')
        tax_df = open_taxids()
        df = pd.read_csv(sys.argv[1], skiprows=1, sep='\t')
        df['seq_rel_date']= pd.to_datetime(df['seq_rel_date'])
        taxid_ls = tax_df['taxid'].to_list()
        abund_ls = tax_df['abund'].to_list()

        for taxid, abund in zip(taxid_ls, abund_ls):
            tmp_df = check_refseq_category(df, taxid)
            if tmp_df.shape[0] == 0:
                tmp_df = check_assembly_level(df, taxid)
            if tmp_df.shape[0] == 0:
                f1.write(f'{taxid}\tna\tna\tna\tna\tna\tna\n')
                print(f'assembly for {taxid} not found')
                continue
            if tmp_df.shape[0] > 1:
                tmp_df = choose_best_assembly(tmp_df)

            best_assembly = tmp_df['ftp_path'].to_list()[0]
            best_assembly = best_assembly + '/' + best_assembly.split('/')[-1] + '_genomic.fna.gz'
            cat = tmp_df['refseq_category'].to_list()[0]
            lev = tmp_df['assembly_level'].to_list()[0]
            rel = tmp_df['release_type'].to_list()[0]
            rep = tmp_df['genome_rep'].to_list()[0]
            date = tmp_df['seq_rel_date'].to_list()[0]
            file_name = os.path.join(sys.argv[3], os.path.basename(best_assembly))
            f1.write(f'{taxid}\t{cat}\t{lev}\t{rel}\t{rep}\t{date}\t{best_assembly}\n')
            f2.write(f'{file_name},{abund}\n')
            subprocess.call(['wget', best_assembly, '-P', sys.argv[3]])


def open_taxids():
    tax_df = pd.read_csv(sys.argv[2], sep='\t')
    tax_df = tax_df[['taxonomy_id', 'fraction_total_reads']]
    tax_df.rename(columns={"taxonomy_id": "taxid", "fraction_total_reads": "abund"},inplace=True)
    return tax_df


def check_refseq_category(df, taxid):
    refseq_ls = ['representative genome', 'reference genome']
    for i in refseq_ls:
        tmp_df = df.loc[((df['taxid'] == taxid) | (df['species_taxid'] == taxid)) & (df['refseq_category'] == i)]
        if tmp_df.shape[0] >= 1:
            break
    return tmp_df


def check_assembly_level(df, taxid):
    assembly_ls = ['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']
    for i in assembly_ls:
        tmp_df = df.loc[((df['taxid'] == taxid) | (df['species_taxid'] == taxid)) & (df['assembly_level'] == i)]
        if tmp_df.shape[0] >= 1:
            break
    return tmp_df

def choose_best_assembly(tmp_df):
    '''
    choose assemblies with:
        major release type
        full genome rep
        newest release date
    '''
    tmp_df = tmp_df.loc[(tmp_df['release_type'] == 'Major') & \
                               (tmp_df['genome_rep'] == 'Full')].copy()
    tmp_df.sort_values(by='seq_rel_date', inplace=True)
    tmp_df = tmp_df[-1:]
    return tmp_df


def download_link():
    pass


if __name__ == '__main__':
    main()
