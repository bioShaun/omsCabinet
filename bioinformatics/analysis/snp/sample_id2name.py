import fire
import pandas as pd


HEADER = ['chrom', 'pos', 'ref', 'alt', 'sample_code']


def id2name(sample_code, sample_order, sample_name):
    sample_code_df = pd.read_csv(
        sample_code, sep='\t', header=None, names=HEADER)
    sample_order_list = [each.strip() for each in open(sample_order)]

    def map_sample_name(sample_code):
        nonlocal sample_order_list
        valid_sample = [each[1] for each in zip(
            sample_code, sample_order_list) if each[0] == '1']
        return ', '.join(valid_sample)

    sample_code_df.loc[:, 'sample_name'] = sample_code_df.sample_code.map(
        map_sample_name)
    sample_code_out_df = sample_code_df.drop('sample_code', axis=1)
    sample_code_out_df.to_csv(sample_name, sep='\t', header=False, index=False)


if __name__ == '__main__':
    fire.Fire(id2name)
