# Covid19Survey

<!-- [![Build Status](https://github.com/andreaskoher/HOPESurvey.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andreaskoher/HOPESurvey.jl/actions/workflows/CI.yml?query=branch%3Amain) -->
<!-- [![Coverage](https://codecov.io/gh/andreaskoher/HOPESurvey.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/andreaskoher/HOPESurvey.jl) -->

This package implements the epidemiological model and data preprocessing steps of the publication **"Monitoring Public behaviour During a Pandemic Using Surveys: Proof-of-Concept Via Epidemic Modelling"** [1]. The epidemiological model is a variant of the popular semi-mechanistic approach of Flaxman et al. [2] and Unwin et. al. [3] as well as the related [Epidemia package](https://imperialcollegelondon.github.io/epidemia/index.html). The code is entirely written in the [Julia language](https://julialang.org/) and we perform inference using the MCMC package [Turing.jl](https://turing.ml/).

# prerequisites

The implementation has been tested on *Julia 1.6*. After starting Julia, please change to the project directory, enter the `pgk>` environment by pressing `]`, activate the project and instantiate to download the correct dependencies:

```
julia> cd("path/to/Covid19Survey")
pkg> activate .
pkg> instantiate
```

# data pre-processing

This package includes the complete pre-processing and analysis pipeline used in the paper *Monitoring Public behaviour During a Pandemic Using Surveys: Proof-of-Concept Via Epidemic Modelling* [1]. If you are not interested in the preprocessing steps and want to perform inference with the same parameters as in the paper, you can just skip this section.

The data folder contains publicly available mobility data from [Google](https://www.google.com/covid19/mobility/), [Apple](https://covid19.apple.com/mobility) and [Danish telco providers](https://covid19.compute.dtu.dk/data-description/telco_data/). It also contains epidemiological data from the [Statens Serum Institute](https://covid19.ssi.dk/) and raw responses to the [HOPE survey](https://hope-project.dk/) on public behaviour and perception of the Covid19 pandemic in Denmark from 01-08-2020 to 01-03-2021.

You can perform the pre-processing steps on the raw moblity data, epidemiological data and behavioural survey data by executing the julia files found in `scripts/mobility`, `scripts/epidata/`, and `scripts/survey`, respectively. We store the pre-processed data as CSV files in the corresponding folders within `data/`. The script `scripts/data_with_predictors.jl` merges mobility, epidemiological and survey data into a final data frame with additional meta-information for the inference model. Similarly, the script `scripts/data_without_predictors.jl` produces a data frame that includes only epidemiological and meta data. Again, we store the final output in the corresponding `data/` folder.

# Inference

We can perform inference conveniently from the terminal by executing either `scripts/inference_with_predictors.jl` or `scripts/inference_without_predictors.jl`:

```
> cd path/to/HOPESurvey/scripts
> julia --threads 5 inference_with_predictors.jl
```

In this case, we start the script with five threads available for parallelization. Both scripts take the following additional parameter:

```
--samples,    -s      number of samples after warm-up (default: 1000)
--chains,     -c      number of chains sampled with multi-threading (default 4)
--warmup,     -w      number of samples to use for warmup/adaptation (default 1000)
--foldername, -f      folder name within HOPESurvey/scripts/ to store results (default \today)
--thinning,   -t      thinning (default 1)
--tree,       -e      maximum tree depth (default 5)
--target,     -r      target acceptance rate (default 0.99)
```

## inference with predictors

We first discuss how to reproduce Fig. 2 and Fig. 3 in the publication using `inference_with_predictors.jl`. The script takes an additional parameter:

```
--predictors, -p    list of predictors from data/data_with_predictors.csv separated by comma (default: total-above18,google,apple,telco)
```

The predictors correspond to columns in the column names in the data file `scripts/data_with_predictors.jl` (see *pre-processing*). As an example, we reproduce the data underlying Fig. 2, i.e., the comparison between self-reported survey data (more precisely: fraction of individuals that report above 18 contacts in total, here: *total-above18*), google mobility data (*google*), apple mobility data (*apple*) and data provided by danish telecommunication companies (*telco*):

```
> cd path/to/HOPESurvey/scripts
> julia --threads 5 inference_with_predictors.jl -s 1000 -c 5 -w 1000 -e 8 -p total-above18,google,apple,telco
```

Similarly, we run the comparison for context depending risk-taking behaviour with a threshold at the 80th percentile (Fig.3). The threshold translates into the following choice of predictors:
```
> cd path/to/HOPESurvey/scripts
> julia --threads 5 inference_with_predictors.jl -s 1000 -c 5 -w 1000 -e 8 -p strangers-above5,family-above3,friends-above4,colleagues-above4
```

## inference without predictors

In Fig. 1 but also in multiple figures in the appendix, we provide a visual comparison between predictors and the inferred reproduction number. We obtain the latter from a daily random walk model of the effective reproduction number without predictors. Executing the following script does the job:
```
> cd path/to/HOPESurvey/scripts
> julia --threads 5 inference_without_predictors.jl -s 1000 -c 5 -w 1000 -e 8
```

[1] A. Koher, F. Jørgensen, M. B. Petersen, S. Lehmann, *Monitoring Public behaviour During a Pandemic Using Surveys: Proof-of-Concept Via Epidemic Modelling*

[2] Flaxman S., Mishra S., Gandy A., Bhatt S. et al. *Estimating the effects of non-pharmaceutical interventions on COVID-19 in Europe*. Nature, 584(7820), 257-261.

[3] Unwin HJT, Mishra S, Bradley VC, Gandy A, et al., *State-level tracking of COVID-19 in the United States.* Nat Commun. 2020 Dec 3;11(1):6189
