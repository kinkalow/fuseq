import os
import shutil

class Checker:
    @staticmethod
    def isdir(d):
        if not os.path.isdir(d):
            print(f'{d} does not exist')
            exit(1)

    @staticmethod
    def isfile(f):
        if not os.path.isfile(f):
            print(f'{f} does not exist')
            exit(1)

    @staticmethod
    def isonefile(li):
        if len(li) != 1:
            print(f'The number of elements in list {li} is {len(li)}, not 1')
            exit(1)

    @staticmethod
    def onshirokane():
        ret = shutil.which('qsub')
        if ret is None:
            print('qsub is not installed')
            exit(1)

    @staticmethod
    def has_blat(path):
        ret = shutil.which(path) if path else None
        if ret is None:
            print('blat is not installed')
            exit(1)
