
# Installation 

**minute** can be installed with conda/mamba via Bioconda channels.

1. Install Conda and Bioconda. Use the
  [Bioconda instructions](https://bioconda.github.io) if you
  don’t have Conda and/or Bioconda, yet.

2. Optionally, but highly recommended is to install [Mamba](https://github.com/mamba-org/mamba), which is a faster alternative to Conda:

        conda install mamba

    If you skip this, write `conda` instead of `mamba` in the next step,
  but the installation will take longer.

3. Install from Bioconda into a new environment:

        mamba create -n minute minute

4. Activate the environment (this must be done before running any of the 
**minute** comands:

        conda activate minute

5. You can check that the installation was successful by running:

        minute --help

This should print **minute** command-line help and exit.