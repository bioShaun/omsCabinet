import fire
import pandas as pd
import pathlib


def main(pcg_cor_file, out_dir=None, cor_cut=0.5, diff_cut=0.3):
    # to purepath object
    pcg_cor_file = pathlib.PurePath(pcg_cor_file)
    out_dir = pathlib.PurePath(out_dir)

    pcg_cor_df = pd.read_table(pcg_cor_file)

    # filter according to cutoff
    f_pcg_cor_df = pcg_cor_df[(pcg_cor_df.module_cor >= cor_cut) &
                              (pcg_cor_df.loc[:, 'diff'] >= diff_cut)]
    f_pcg_cor_file = pcg_cor_file.with_suffix('.filter.txt')
    f_pcg_cor_df.to_csv(f_pcg_cor_file, sep='\t', index=False)

    # output module gene list
    gene_type = pcg_cor_file.name.split('.')[0]
    modules = f_pcg_cor_df.module_name.unique()
    for each_module in modules:
        each_module_df = f_pcg_cor_df[f_pcg_cor_df.module_name == each_module]
        module_list_df = each_module_df.loc[:, 'mRNA']
        each_module_file = out_dir / f'{gene_type}.{each_module}.gene'
        module_list_df.to_csv(each_module_file, sep='\t',
                              index=False, header=False)

    # module mRNA count
    mc_df = pd.DataFrame(f_pcg_cor_df.groupby(['module_name']).size())
    mc_df.columns = ['co_exp_pcg']
    mc_file = out_dir / f'{gene_type}.module.co_mRNA.num'
    mc_df.to_csv(mc_file, sep='\t')


if __name__ == '__main__':
    fire.Fire(main)
