import os
import shutil
from setuptools import setup
from fuseq import __version__

scripts_dir = 'build/_scripts'
if not os.path.exists(scripts_dir):
    os.makedirs(scripts_dir)
shutil.copyfile('fuseq.py', f'{scripts_dir}/fuseq')

setup(
    name='fuseq',
    version=__version__,
    description='Fusion gene sequences are obtained',
    packages=['fuseq'],
    scripts=[f'{scripts_dir}/fuseq'],
    python_requires='>=3.6',
)
