#! /usr/bin/env python


import fire
import requests
from bs4 import BeautifulSoup
import urllib.parse
import pathlib
import asyncio
import aiohttp
import time


URL = 'http://cantata.amu.edu.pl/download.php'
CURRENT_DIR = pathlib.Path().cwd()

SEMA = asyncio.Semaphore(5)


async def get_file(url, res_dict):
    with (await SEMA):
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as resp:
                file_name = pathlib.PurePath(url).name
                print('Downloading {n}...'.format(n=file_name))
                url_content = await resp.text()
                res_dict[url] = url_content
                time.sleep(10)


def download_files(out_dir=CURRENT_DIR,
                   url=URL):
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'lxml')
    download_hrefs = soup.find(
        'table', {'class': 'table table-hover'}).findAll('a')
    href_dict = dict()
    for each_href in download_hrefs:
        if 'href' not in each_href.attrs:
            continue
        each_url = urllib.parse.urljoin(url,
                                        each_href['href'])
        each_name = pathlib.PurePath(each_url).name
        each_file = pathlib.PurePath(out_dir) / each_name
        if pathlib.Path(each_file).exists():
            print('File {n} exists, skip download.'.format(
                n=each_name
            ))
        else:
            href_dict[each_url] = each_file
    res_dict = dict()
    if href_dict:
        tasks = [get_file(host, res_dict) for host in href_dict]
        # # for test
        # tasks = tasks[:8]
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.wait(tasks))
        for each_url in res_dict:
            each_file = href_dict[each_url]
            each_content = res_dict[each_url]
            with open(each_file, 'w') as file_inf:
                file_inf.write(each_content)
    else:
        print('All file is downloaded.')


if __name__ == '__main__':
    fire.Fire(download_files)
