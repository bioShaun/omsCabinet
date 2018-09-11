import itertools
import sys
import asyncio
import aioftp
import pathlib
import fire
import re
import pandas as pd
import gzip
import concurrent


SEMA = asyncio.Semaphore(5)
CURRENT_DIR = pathlib.Path().cwd()
HOST = "ftp.ncbi.nlm.nih.gov"


async def get_inf(host, name=None, retry_lim=1):
    file_list = []
    retry_num = 0
    with (await SEMA):
        while retry_num <= retry_lim:
            try:
                async with aioftp.ClientSession(
                        host, socket_timeout=10) as client:
                    path = '/blast/db/'
                    for path, info in (await client.list(path)):
                        if info['type'] == 'file' and (path.suffix == '.gz' or
                                                       path.suffix == '.md5'):
                            if name is None:
                                name = '.*'
                            pattern = re.compile(
                                '{pref}.(\w+).tar.gz.*'.format(
                                    pref=name
                                ))
                            if pattern.match(path.name):
                                file_list.append(path)
                    return file_list
            except concurrent.futures._base.TimeoutError:
                retry_num += 1
        return None


async def get_file(host, path, out_dir=CURRENT_DIR, retry_lim=2):
    with (await SEMA):
        retry_times = 0
        while retry_times <= retry_lim:
            try:
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
                    return True
            except TimeoutError:
                retry_times += 1
        return False


def download_ncbi_blastdb(database,
                          out_dir=CURRENT_DIR,
                          test=False):

    loop = asyncio.get_event_loop()
    inf_tasks = [
        get_inf(HOST, database)
    ]
    file_list = loop.run_until_complete(asyncio.wait(inf_tasks))
    print(file_list)
    print(file_list.result)
    loop.close()
    # download_tasks = [
    #     get_file(HOST, each_path) for each_path
    #     in file_list]
    # download_status = loop.run_until_complete(asyncio.wait(download_tasks))
    # loop.close()

    # total_works = len(download_status)
    # success_works = sum(download_status)
    # failed_works = total_works - success_works
    # print(f'{total_works} files to be downloaded.')
    # print(f'{success_works} success.')
    # print(f'{failed_works} failed.')


@asyncio.coroutine
def spin(msg):
    write, flush = sys.stdout.write, sys.stdout.flush
    for char in itertools.cycle('|/-\\'):
        status = char + ' ' + msg
        write(status)
        flush()
        write('\x08' * len(status))
        try:
            yield from asyncio.sleep(.1)
        except asyncio.CancelledError:
            break
    write(' ' * len(status) + '\x08' * len(status))


def list_db():

    @asyncio.coroutine
    def get_db_inf(msg='Fetching FTP information!'):
        spinner = asyncio.async(spin(msg))
        file_list = yield from get_inf(HOST)
        spinner.cancel()
        return file_list

    loop = asyncio.get_event_loop()
    file_list = loop.run_until_complete(get_db_inf())
    loop.close()
    if file_list is None:
        print('Failed to connect to the FTP host.')
        print(f'Please check host IP [{HOST}] and try again!')
        sys.exit(1)
    file_set = {each_path.stem.split('.')[0] for
                each_path in file_list}
    file_list = sorted(list(file_set))
    db_title = 'NCBI Blast database'
    print('=' * 30)
    print(f'|{db_title:^28}|')
    print('=' * 30)
    for each_db in file_list:
        print(f'|{each_db:^28}|')
    print('=' * 30)


if __name__ == '__main__':
    fire.Fire()
