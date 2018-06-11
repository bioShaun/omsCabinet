#! /usr/bin/python3
# coding=utf-8


# change all files and directorys from group a TO b
# execute this script with root account

import os
import pathlib
import fire

BROKEN_LNK = "链接失效，文件不存在"


def check_broken_symlink(path):
    path = pathlib.Path(path)
    if path.is_symlink() and (not path.exists()):
        print("{bl_chn}: {pt}".format(
            bl_chn=BROKEN_LNK,
            pt=path))
        return True
    else:
        return False


def change_files_group(dir_path, old_group='sycai', new_group='xms'):
    for root, dirs, files in os.walk(dir_path):
        dir_file_list = [dir_path]
        dir_file_list.extend(dirs)
        dir_file_list.extend(files)
        for each_path in dir_file_list:
            each_path = pathlib.PurePath(root) / each_path
            if not check_broken_symlink(each_path):
                if pathlib.Path(each_path).group() == old_group:
                    os.system('chgrp {ng} {pt}'.format(ng=new_group,
                                                       pt=each_path))


if __name__ == '__main__':
    fire.Fire(change_files_group)
