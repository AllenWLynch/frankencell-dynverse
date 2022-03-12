#!/bin/bash

snakemake -p -s frankencell/pipeline/SnakeFile \
	--configfiles frankencell/pipeline/2022-03-08_experiment/2022-03-08_experiment.yaml frankencell/pipeline/2022-03-08_experiment/2022-03-08_experiment_no_MIRA.yaml \
	--config conda_env=benchmarking_python restart=False \
	--use-conda \
	--use-envmodules \
	-j 50 \
	--resources mem_mb=100000 cpus=50 \
	--profile ~/.config/snakemake/allen
