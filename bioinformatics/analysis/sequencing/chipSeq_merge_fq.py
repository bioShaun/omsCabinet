import re
import fire
import asyncio
from pathlib import Path

clean_reads_pattern = '(.*).R([1,2]).clean.fastq.gz'
low_qual_reads_pattern = '(.*).R([1,2]).lowqual.fastq.gz'
adapter_reads_pattern = '(.*).R([1,2]).adapter.fastq.gz'


async def run(cmd, semaphore):
    with await semaphore:
        proc = await asyncio.create_subprocess_shell(
            cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)
        stdout, stderr = await proc.communicate()
        # await asyncio.sleep(3)
        # print(stdout)
    return stdout, stderr


def load_fq_files(chip_fq_dir, reads_pattern, fq_dict):
    for each_file in chip_fq_dir.iterdir():
        re_test = reads_pattern.match(each_file.name)
        if re_test:
            sample_name, reads_num = re_test.groups()
            try:
                fq_dict[sample_name][reads_num].append(str(each_file))
            except KeyError:
                fq_dict.setdefault(sample_name,
                                   {})[reads_num] = [str(each_file)]


def merge_chip_fq(chip_fq_dir, thread=4):
    fq_dict = dict()
    chip_fq_dir = Path(chip_fq_dir)
    clean_fq_re = re.compile(clean_reads_pattern)
    load_fq_files(chip_fq_dir, clean_fq_re, fq_dict)
    other_fq_dir = chip_fq_dir / '002'
    low_qual_fq_re = re.compile(low_qual_reads_pattern)
    load_fq_files(other_fq_dir, low_qual_fq_re, fq_dict)
    adapter_fq_re = re.compile(adapter_reads_pattern)
    load_fq_files(other_fq_dir, adapter_fq_re, fq_dict)
    merge_dir = chip_fq_dir / 'rawdata'
    merge_dir.mkdir(exist_ok=True)
    cmd_list = list()
    md5_cmd_list = list()
    for each_sample in fq_dict:
        for each_read in fq_dict[each_sample]:
            reads_files = ' '.join(fq_dict[each_sample][each_read])
            merged_read = merge_dir / f'{each_sample}.R{each_read}.fq.gz'
            md5_file = merge_dir / f'{each_sample}.R{each_read}.fq.gz.md5'
            if not merged_read.exists():
                merge_cmd = f'cat {reads_files} > {merged_read}'
                cmd_list.append(merge_cmd)
            if not md5_file.exists():
                md5_cmd = f'md5sum {merged_read} > {md5_file}'
                md5_cmd_list.append(md5_cmd)
    semaphore = asyncio.Semaphore(thread)
    loop = asyncio.get_event_loop()
    if cmd_list:
        coro_list = [run(cmd, semaphore) for cmd in cmd_list]
        loop.run_until_complete(asyncio.wait(coro_list))
    if md5_cmd_list:
        md5_coro_list = [run(cmd, semaphore) for cmd in md5_cmd_list]
        loop.run_until_complete(asyncio.wait(md5_coro_list))
    loop.close()


if __name__ == '__main__':
    fire.Fire(merge_chip_fq)
