from setuptools import setup, find_packages

setup(name='aavolve',
      version='0.1.0',
      description='AAVovle is a tools for analysing sequencing data from directed evolution experiments with shuffled libraries.',
      author='Suzanne Scott',
      url='https://github.com/szsctt/AAVolve',
      packages=find_packages(),
      install_requires=['biopython', 'pysam', 'numpy', 'scipy', 'pandas', 'tqdm'])