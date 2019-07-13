import click
from pathlib import Path
from omutils import async_batch_sh_jobs
from collections import defaultdict

CWD = Path().cwd()


@click.command()
@click.option('--fq_dir', '-f', multiple=True)
@click.option('--out_dir', '-o', default=CWD, type=click.Path())
@click.option('--fq_suffix', '-s', default='gz', type=click.STRING)
@click.option('--mode',
              '-m',
              default='ln',
              type=click.Choice(['ln', 'cp']),
              help='ln or cp fq files not need to combine.')
@click.option('--thread', '-t', default=4, type=click.INT)
def main(fq_dir, out_dir, fq_suffix, mode, thread):
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    fq_dict = defaultdict(list)
    for fq_dir_i in fq_dir:
        fq_files = Path(fq_dir_i).glob(f'*{fq_suffix}')
        for fq_i in fq_files:
            fq_dict[fq_i.name].append(fq_i)

    cmd_list = []
    for fq_i in fq_dict:
        fq_i_paths = fq_dict[fq_i]
        if len(fq_i_paths) == 1:
            if mode == 'ln':
                mode = 'ln -s'
            cmd = f'{mode} {fq_i_paths[0]} {out_dir}/{fq_i}'
        else:
            fq_i_paths = [str(path_i) for path_i in fq_i_paths]
            fq_i_path_str = ' '.join(fq_i_paths)
            cmd = f'cat {fq_i_path_str} > {out_dir}/{fq_i}'
        cmd_list.append(cmd)

    async_batch_sh_jobs(cmd_list, thread)


if __name__ == "__main__":
    main()
