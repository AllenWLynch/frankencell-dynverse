# Frankencell installation and usage

This repository contains the code needed to reproduce the data used in benchmarking [MIRA](https://www.biorxiv.org/content/10.1101/2021.12.06.471401v1.full.pdf).

### Installation

First, set up a directory:

```
$ mkdir ~/frankencell 
$ cd ~/frankencell
```

Next, download the Frankencell source code from github and set permissions to use the "franken-cells" script.

```
$ git clone https://github.com/AllenWLynch/frankencell-dynverse.git
$ cd frankencell-dynverse
$ chmod u+x franken-cells
```

Next, install the required packages from conda and activate the environment:

```
$ conda create -n mira-benchmarking -c conda-forge -c bioconda numpy scipy pandas anndata networkx joblib tqdm snakemake-minimal
$ conda activate mira-benchmarking
```

And also install the "dynclipy" package from github using PIP:

```
$ pip install git+https://github.com/dynverse/dynclipy.git
```

Finally, set up your R installation with the required packages from dynverse:

```
$ R
> install.packages("devtools")
> devtools::install_github("dynverse/dyno")
> devtools::install_github("dynverse/dynwrap")
> devtools::install_github("dynverse/dynutils")
```

### Usage

To run the pipeline, make the following directorites:

```
$ mkdir data
$ mkdir data/mira_benchmarking
```

Then go to [Zenodo](https://zenodo.org/record/6390740#.YkHgaO7MIUF) and download "startdata_rna.h5ad", "startdata_atac.h5ad", and "data_config.yaml". Place these three files in your working directory alongside the "franken-cells" script.

The command "franken-cells gen" starts a snakemake pipeline. Pass the path to "data_config.yaml" as the first argument, then give your available memory in megabytes (this determines how many of certain types of jobs may run concurrently).
After the "--snake-args" flag, every additional argument is passed directly to the [snakemake pipeline command-line executor](https://snakemake.readthedocs.io/en/stable/executing/cli.html). The command below runs a dry run of the 
pipeline to make sure everything is set up correctly.

```
$ ./franken-cells gen data_config.yaml --mem-mb {your available mem here in MB} --restart False --snake-args --keep-going --cores {your available cores here} --rerun-incomplete -p -n
```

To execute the pipeline, remove the "-n" flag from the end of the above command.
