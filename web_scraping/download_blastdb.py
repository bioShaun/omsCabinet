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
import logging
from tabulate import tabulate
from collections import Counter


def init_logger(log_name=None, log_file=None, level=logging.INFO,
                streamhandler=True):
    if log_name is None:
        logger = logging.getLogger(__name__)
    else:
        logger = logging.getLogger(log_name)
    logger.setLevel(level=level)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # FileHandler
    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    # StreamHandler
    if streamhandler:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)
    return logger


DEFAULT_LOGGER = init_logger()


SEMA = asyncio.Semaphore(5)
CURRENT_DIR = pathlib.Path().cwd()
HOST = "ftp.ncbi.nlm.nih.gov"


async def get_inf(host, name=None,
                  retry_lim=1,
                  logger=DEFAULT_LOGGER):
    file_list = []
    retry_num = 0
    out_name = name
    # name is None means to retrieve all files of blastdb
    if name is None:
        out_name = 'NCBI Blast Database'
    with (await SEMA):
        while retry_num <= retry_lim:
            try:
                async with aioftp.ClientSession(
                        host, socket_timeout=10) as client:
                    path = '/blast/db/'
                    logger.info(
                        f'Retriving blastdb files. Try {retry_num + 1}.')
                    for path, info in (await client.list(path)):
                        if info['type'] == 'file' and (path.suffix == '.gz' or
                                                       path.suffix == '.md5'):
                            if name is None:
                                name = '.*'
                            pattern = re.compile(
                                '{pref}\..*tar.gz.*'.format(
                                    pref=name
                                ))
                            if pattern.match(path.name):
                                file_list.append(path)
                    file_number = len(file_list)
                    logger.info(
                        f'Total {file_number} files for [{out_name}] database.')
                    return file_list
            except concurrent.futures._base.TimeoutError:
                retry_num += 1
            except aioftp.errors.StatusCodeError:
                retry_num += 1
            except ConnectionResetError:
                retry_num += 1
        logger.error('Failed to connect to the FTP host.')
        logger.error(f'Please check host IP [{host}] and try again!')
        return None


async def get_file(host, path, out_dir=CURRENT_DIR,
                   retry_lim=2, logger=DEFAULT_LOGGER):
    with (await SEMA):
        retry_times = 0
        while retry_times <= retry_lim:
            try:
                async with aioftp.ClientSession(
                        host, socket_timeout=30) as client:
                    # get download file stat
                    if await client.exists(path):
                        stat = await client.stat(path)
                        size = int(stat["size"])
                    else:
                        logger.error(f'{file_name} not exists in server!')
                    file_name = pathlib.PurePath(path).name
                    outfile = out_dir / file_name
                    if pathlib.Path(outfile).exists():
                        outfile_stat = pathlib.Path(outfile).stat()
                        outfile_size = outfile_stat.st_size
                    else:
                        outfile_size = 0
                    if outfile_size != size:
                        logger.info(f'|Downloading...| {file_name}')
                        file_out = open(outfile, 'ab')
                        async with client.download_stream(
                                path, offset=outfile_size) as stream:
                            async for block in stream.iter_by_block():
                                file_out.write(block)
                        stream.close()
                        file_out.close()
                    else:
                        pass
                    logger.info(f'|Downloaded| {file_name}')
                    return True
            except TimeoutError:
                retry_times += 1
            except concurrent.futures._base.TimeoutError:
                retry_times += 1
            except aioftp.errors.StatusCodeError:
                retry_times += 1
            except ConnectionResetError:
                retry_times += 1
        logger.error(f'|Failed!| {file_name}')
        return False


@asyncio.coroutine
def spin(msg):
    '''
    show a spinning line
    '''
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


@asyncio.coroutine
def get_db_inf(msg='Fetching FTP information!',
               show_spin=True,
               db=None):
    # when downloading, don not show spinning line
    if show_spin:
        spinner = asyncio.async(spin(msg))
        file_list = yield from get_inf(HOST, name=db)
        spinner.cancel()
    else:
        file_list = yield from get_inf(HOST, name=db)
    return file_list


def file_counter(file_list, file_suffix, file_label):
    db_file_list = [each_path.stem.split('.')[0] for
                    each_path in file_list
                    if each_path.suffix == file_suffix]
    db_file_count = Counter(db_file_list)
    db_file_dict = {}
    db_file_dict['NCBI_Blast_Database'] = list(db_file_count.keys())
    db_file_dict[file_label] = list(db_file_count.values())
    db_file_df = pd.DataFrame(db_file_dict)
    return db_file_df.loc[:, ['NCBI_Blast_Database', file_label]]


def download_ncbi_blastdb(database,
                          out_dir=CURRENT_DIR,
                          test=False):

    logger_file = pathlib.PurePath(out_dir) / 'download.log.txt'
    dl_logger = init_logger(log_name='download_ncbi_blastdb',
                            log_file=logger_file)

    loop = asyncio.get_event_loop()
    db_msg = f'Fetching dababase [{database}] files.'
    file_list = loop.run_until_complete(
        get_db_inf(db_msg,
                   db=database,
                   show_spin=False))

    if file_list is None:
        return

    if test:
        for each_file in file_list:
            print(each_file)
        loop.close()
    else:
        download_tasks = [
            get_file(HOST, path=each_path, logger=dl_logger) for each_path
            in file_list]
        download_status, _ = loop.run_until_complete(
            asyncio.wait(download_tasks))
        loop.close()
        total_works = len(download_status)
        success_works = sum([each.result() for each in download_status])
        failed_works = total_works - success_works
        dl_logger.info(f'{total_works} files to be downloaded.')
        dl_logger.info(f'{success_works} success.')
        dl_logger.info(f'{failed_works} failed.')
        if failed_works:
            dl_logger.info('Check log file for failed files.')


def list_db():
    loop = asyncio.get_event_loop()
    file_list = loop.run_until_complete(get_db_inf())
    loop.close()
    if file_list is None:
        return
    db_file_list = [each_path.stem.split('.')[0] for
                    each_path in file_list
                    if each_path.suffix == '.gz']
    db_file_num = file_counter(file_list, '.gz', 'DB_File_Number')
    db_md5_num = file_counter(file_list, '.gz', 'DB_MD5_Number')
    total_db_files = db_file_num.merge(db_md5_num)
    total_db_files = total_db_files.sort_values(['NCBI_Blast_Database'])
    total_db_files.index = list(range(len(total_db_files)))
    print(tabulate(total_db_files, headers='keys', tablefmt='psql'))


if __name__ == '__main__':
    fire.Fire({
        'list_db': list_db,
        'download_ncbi_blastdb': download_ncbi_blastdb,
    })
