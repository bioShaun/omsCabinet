import asyncio
import aioftp
import pathlib
import fire
import re
import pandas as pd
import gzip


SEMA = asyncio.Semaphore(5)
CURRENT_DIR = pathlib.Path().cwd()
HOST = "ftp.ensemblgenomes.org"


def ensembl_file_path(species, g_version, d_version):
    base_path = pathlib.PurePath('/pub/plants/release-{ver}/'.format(
        ver=d_version
    ))
    file_pref = '{sp}.{gv}'.format(
        sp=species.capitalize(), gv=g_version
    )
    cds_file_name = '{pref}.cds.all.fa.gz'.format(
        pref=file_pref
    )
    cds_file = base_path / 'fasta/{sp}/cds/{name}'.format(
        sp=species, name=cds_file_name
    )
    gtf_file_name = '{pref}.{dv}.gtf.gz'.format(
        pref=file_pref, dv=d_version
    )
    gtf_file = base_path / 'gtf/{sp}/{name}'.format(
        sp=species, name=gtf_file_name
    )
    pep_file_name = '{pref}.pep.all.fa.gz'.format(
        pref=file_pref
    )
    pep_file = base_path / 'fasta/{sp}/pep/{name}'.format(
        sp=species, name=pep_file_name
    )
    test_file_name = '{pref}.dna.toplevel.fa.gz.fai'.format(pref=file_pref)
    test_file = base_path / 'fasta/{sp}/dna_index/{name}'.format(
        sp=species, name=test_file_name
    )
    out_dict = {
        'cds': cds_file,
        'gtf': gtf_file,
        'pep': pep_file,
        'test_file': test_file
    }
    return out_dict


async def get_inf(host, species, version, info_dict):
    with (await SEMA):
        async with aioftp.ClientSession(host) as client:
            path = '/pub/plants/release-{ver}/fasta/{sp}/cds/'.format(
                ver=version, sp=species
            )
            for path, info in (await client.list(path)):
                if info['type'] == 'file' and path.suffix == '.gz':
                    pattern = re.compile('{pref}.(\S+).cds.all.fa.gz'.format(
                        pref=species.capitalize()
                    ))
                    g_version = pattern.match(path.name).groups()[0]
                    info_dict[species] = g_version


async def get_file(host, path, out_dir=CURRENT_DIR):
    with (await SEMA):
        async with aioftp.ClientSession(host, socket_timeout=30) as client:
            # get download file stat
            if await client.exists(path):
                stat = await client.stat(path)
                size = int(stat["size"])
            else:
                print('{fi} not exists!'.format(fi=path))
            file_name = pathlib.PurePath(path).name
            outfile = out_dir / file_name
            if pathlib.Path(outfile).exists():
                outfile_stat = pathlib.Path(outfile).stat()
                outfile_size = outfile_stat.st_size
                if outfile_size != size:
                    print('Downloading {fi}'.format(fi=file_name))
                    file_out = gzip.open(outfile, 'ab')
                    async with client.download_stream(
                            path, offset=outfile_size) as stream:
                        file_out.write(stream)
                    file_out.close()
                else:
                    pass
            else:
                print('Downloading {fi}'.format(fi=file_name))
                await client.download(path, outfile, write_into=True)
            print('Download {fi} finished.'.format(fi=file_name))


def download_ensembl_files(species_file, version,
                           out_dir=CURRENT_DIR, test=False):
    species_df = pd.read_table(species_file, header=None,
                               names=['species'])

    def format_sp(sp):
        return '_'.join(sp.split()).lower()

    species_df.loc[:, 'species'] = species_df.species.map(format_sp)

    info_dict = dict()
    loop = asyncio.get_event_loop()
    inf_tasks = [
        get_inf(HOST, each_sp, version, info_dict) for
        each_sp in species_df.species]
    loop.run_until_complete(asyncio.wait(inf_tasks))

    # file_type = ['cds', 'gtf']
    file_type = ['pep']
    if test:
        file_type = ['test_file']

    file_list = list()
    for each_sp in species_df.species:
        each_sp_gv = info_dict[each_sp]
        each_sp_files = ensembl_file_path(each_sp, each_sp_gv,
                                          version)
        for each_file in file_type:
            file_list.append(each_sp_files[each_file])

    download_tasks = [
        get_file(HOST, each_path) for each_path
        in file_list]
    loop.run_until_complete(asyncio.wait(download_tasks))
    loop.close()


if __name__ == '__main__':
    fire.Fire(download_ensembl_files)
