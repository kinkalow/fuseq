import os
import subprocess

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
    def isshirokane():
        cmd = 'type qsub > /dev/null'
        p = subprocess.Popen(cmd, shell=True)
        p.communicate()
        if p.returncode:
            print('Not in Shirokane')
            exit(1)

    @staticmethod
    def has_tools():
        cmd = 'type blat > /dev/null'
        p = subprocess.Popen(cmd, shell=True)
        p.communicate()
        if p.returncode:
            print('blat is not installed')
            exit(1)
