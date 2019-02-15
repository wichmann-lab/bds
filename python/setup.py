from distutils.core import setup
setup(
    name="bds",
    version="0.1",
    packages=['bds'],
    data_files = [('share/models', ['../models/bds.stan', '../models/bds_lps.stan'])],
    author='Bernhard Lang',
    author_email='bernhard.lang@student.uni-tuebingen.de',
    url='https://www.github.com/wichmann-lab/bds'
)
