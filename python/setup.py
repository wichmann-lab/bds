from setuptools import setup, find_packages
setup(
    name="bds",
    version="0.2",
    packages=find_packages(),
    install_requires=['pystan', 'pandas', 'dfply', 'plotnine'],
    author='Bernhard Lang',
    author_email='bernhard.lang@student.uni-tuebingen.de',
    url='https://www.github.com/wichmann-lab/bds'
)
