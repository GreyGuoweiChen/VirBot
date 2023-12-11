from setuptools import setup, find_packages
import glob
import os
import pkg_resources
from setuptools.command.install import install as _install
import zipfile

from virbot import __version__

class install(_install):
    def pull_first(self):
        """This script is in a git directory that can be pulled."""
        gitdir = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(gitdir, "virbot", "data", "ref")
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        if os.listdir(data_dir):
            return

        zip_data = os.path.join(gitdir, "virbot", "data", "ref.zip")
        if not os.path.exists(zip_data):
            import git
            os.chdir(gitdir)
            g = git.cmd.Git(gitdir)
            try:
                g.execute(['git', 'lfs', 'install'])
                g.execute(['git', 'lfs', 'fetch'])
                g.execute(['git', 'lfs', 'pull'])
            except git.exc.GitCommandError:
                print("Warning git-lfs is not installed - please manually download and unzip the reference files!")
            os.chdir(cwd)

        with zipfile.ZipFile(zip_data,"r") as zip_ref:
            zip_ref.extractall() 

    def run(self):
        self.pull_first()
        super().run()

setup(name='virbot',
      version=__version__,
      packages=find_packages(),
      package_data={"virbot":["data/ref/*"]},
      install_requires=[],
      setup_requires=[],
      description='VirBot: RNA viral contig detector for metagenomic data',
      url='https://github.com/GreyGuoweiChen/RNA_virus_detector.git',
      author='Guowei Chen',
      author_email='',
      entry_points={"console_scripts": ["virbot = virbot.VirBot:main"]},
      include_package_data=True,
      keywords=[],
      cmdclass={'install': install},
      zip_safe=False)
