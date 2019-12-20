import fire
import pandas as pd


def format_df(table_file, by, outfile, sep='||', empty_rep=None):
    '''
    condense multi row to one by identical column value
    '''

    gene_df = pd.read_table(table_file)
    gene_df.dropna(inplace=True)

    def my_unique(x, sep='||'):
        nonlocal empty_rep
        unique_x = pd.unique(x.dropna())
        if str(unique_x.dtype) == 'float64':
            unique_x = unique_x.astype('int')
        unique_x = [str(each) for each in unique_x]
        if not unique_x:
            if empty_rep is None:
                return None
            else:
                unique_x = [empty_rep]
        try:
            sep.join(unique_x)
        except TypeError:
            print(unique_x)
            print(x)
        return sep.join(unique_x)

    gene_df = gene_df.groupby(by).agg(my_unique)
    gene_df.to_csv(outfile, sep='\t')


if __name__ == '__main__':
    fire.Fire(format_df)
