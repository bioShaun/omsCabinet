import fire
import pandas as pd


HEADER = [
    'uniprot_id',
    'ensembl_id',
    'interpro_id',
    'interpro_description',
]


def extract_uniprot_des(uniprot_text, outfile):
    uniprot_des_dict = dict()
    ensembl_plants_list = list()
    interpro_id_list = list()
    interpro_des_list = list()

    with open(uniprot_text) as uniprot_text_inf:
        for eachline in uniprot_text_inf:
            if eachline.startswith('AC'):
                each_id = eachline.split()[-1]
                each_id = each_id.split(';')[0].strip()
            elif eachline.startswith('DR') and 'EnsemblPlants;' in eachline:
                ensembl_id = eachline.split(';')[-1].strip().rstrip('.')
                ensembl_plants_list.append(ensembl_id)
            elif eachline.startswith('DR') and 'InterPro;' in eachline:
                interpro_id_list.append(eachline.split(';')[1].strip())
                interpro_des_list.append(
                    eachline.split(';')[2].strip().rstrip('.'))
            elif eachline.startswith('//'):
                ensembl_id = '|'.join(ensembl_plants_list)
                interpro_id = '|'.join(interpro_id_list)
                interpro_description = '|'.join(interpro_des_list)
                uniprot_des_dict.setdefault('uniprot_id', []).append(each_id)
                uniprot_des_dict.setdefault(
                    'ensembl_id', []).append(ensembl_id)
                uniprot_des_dict.setdefault(
                    'interpro_id', []).append(interpro_id)
                uniprot_des_dict.setdefault(
                    'interpro_description', []).append(interpro_description)
                ensembl_plants_list = list()
                interpro_id_list = list()
                interpro_des_list = list()
            else:
                pass
    uniprot_des_df = pd.DataFrame(uniprot_des_dict)
    uniprot_des_df.to_csv(outfile, sep='\t', index=False, columns=HEADER,
                          na_rep='--')


if __name__ == '__main__':
    fire.Fire(extract_uniprot_des)
