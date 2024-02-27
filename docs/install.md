
# Installation 

**minute** can be installed with conda via Bioconda channels following these
steps:

1. If you don't have Conda and Bioconda: install them following the
  [Bioconda instructions](https://bioconda.github.io).

2. Create a new environment from Bioconda:

        conda create -n minute minute

3. Activate the environment (this must be done before running any of the 
**minute** commands):

        conda activate minute

4. You can check that the installation was successful by running:

        minute --help

This should print **minute** command-line help and exit.


If you are a developer: Instead of the above, install the dependencies with
`conda env create -n minute -f environment.yml` and then install **minute**:
`pip install .`.
