import re
import attr
import click
import asyncio
from pathlib import Path
from collections import defaultdict

READS_PAT = re.compile(r'(\S+.R[1,2]).*fastq.gz')


@attr.s(hash=True)
class readsInf:
    read_name = attr.ib()
    read_stat = attr.ib()
    read_path = attr.ib(hash=False, cmp=False)


async def async_sh_job(cmd, sema):
    with (await sema):
        p = await asyncio.create_subprocess_shell(
            cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT)
        return (await p.communicate())[0].splitlines()


def async_batch_sh_jobs(cmd_list, thread=2):
    if cmd_list:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        loop = asyncio.get_event_loop()
        semaphore = asyncio.Semaphore(thread)
        coro_list = [async_sh_job(cmd, semaphore) for cmd in cmd_list]
        try:
            loop.run_until_complete(asyncio.wait(coro_list))
        finally:
            loop.close()


def fq_file_cls(fq_dir, fq_dict):
    fq_dir = Path(fq_dir).resolve()
    for file_i in fq_dir.glob('**/*'):
        if '$RECYCLE.BIN' in file_i.parent.parts:
            continue
        read_i = READS_PAT.match(file_i.name)
        if read_i:
            read_i_name = read_i.groups()[0]
            read_i_obj = readsInf(read_i_name, file_i.stat().st_size, file_i)
            if read_i_name in fq_dict:
                if read_i_obj in fq_dict[read_i_name]:
                    pass
            fq_dict[read_i_name].append(read_i_obj)


@click.command()
@click.option('-f',
              '--fq_dir',
              help='fq file dir.',
              type=click.Path(exists=True, file_okay=False),
              multiple=True,
              required=True)
@click.option('-o',
              '--out_dir',
              help='fq out dir.',
              type=click.Path(),
              required=True)
@click.option('-t', '--thread', help='thread use.', type=click.INT, default=1)
def main(fq_dir, out_dir, thread):
    out_dir = Path(out_dir).resolve()
    log_file = out_dir / 'log.txt'
    log_inf = log_file.open('a')
    out_dir.mkdir(parents=True, exist_ok=True)
    fq_dict = defaultdict(list)
    cmd_list = []
    for dir_i in fq_dir:
        fq_file_cls(dir_i, fq_dict)
    for read_i in fq_dict:
        read_i_inf = fq_dict[read_i]
        read_i_out = out_dir / f'{read_i}.raw.fastq.gz'
        if len(read_i_inf) == 1:
            cmd = f'ln -s {read_i_inf[0].read_path} {read_i_out}'
        else:
            reads_p_list = sorted(
                [str(read_i.read_path) for read_i in read_i_inf])
            reads_p_str = ' '.join(reads_p_list)
            cmd = f'cat {reads_p_str} > {read_i_out}'
        cmd_list.append(cmd)
        print(cmd, file=log_inf)
    log_inf.close()
    async_batch_sh_jobs(cmd_list, thread)


if __name__ == "__main__":
    main()
