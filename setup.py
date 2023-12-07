from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from virbot import __version__

setup(name='virbot',
      version=__version__,
      packages=find_packages(),
      package_data={"virbot":["data/ref/*"]},
      install_requires=[
            'pandas>=1.0.1'
        ],
      description='VirBot: RNA viral contig detector for metagenomic data',
      url='https://github.com/GreyGuoweiChen/RNA_virus_detector.git',
      author='Guowei Chen',
      author_email='',
      entry_points={"console_scripts": ["virbot = virbot.VirBot:main"]},
      include_package_data=True,
      keywords=[],
      zip_safe=False)
