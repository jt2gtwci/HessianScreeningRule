<!-- README.md is generated from README.Rmd. Please edit that file -->

# Code for the Hessian Screening Rule

<!-- badges: start -->
<!-- badges: end -->

## Results

The results from the simulations are stored in the [results
folder](results/). The figures and tables in the paper, generated from
these results, are stored in [`figures/`](figures/) and
[`tables/`](tables/) respectively.

## Reproducing the Results

To reproduce the results, we recommend you use the singularity
container. See the release section on GitHub and download the container
from there. To run an experiment from the singularity container, call

    singularity run --no-home --bind results:/project/results container.sif <script>

where `<script>` should be a name of a script to run from the
[experiments folder](experiments/) folder, such as
`experiments/simulateddata.R`. The results will then be output to the
`results` folder.

### Re-building the Singularity Container

If you want to re-build the singularity container (or simply want to
clone the repo to your local drive), you can do so via the following
steps.

1.  Clone the repository to your local hard drive. On linux, using SSH
    authentication, run

        git clone git@github.com:jt2gtwci/HessianScreeningRule.git

2.  Navigate to the root of the repo and build the singularity container
    by calling

        cd HessianScreeningRule
        sudo singularity build container.sif Singularity

Then proceed as in [Reproducing the Results](#reproducing-the-results)
to run the experiments.

### Running Experiments without Singularity

Alternatively, you may also reproduce the results by cloning this
repository, then either opening the `HessianScreening.Rproj` file in R
Studio or starting R in the root directory of this folder (which will
activate the renv repository) and then run

    renv::restore()

to restore the project library. Then build the R package (see below) and
run the simulations directly by running the scripts in the experiments
folder. This is not recommended, however, since it, unlike the
Singularity container approach, does not exactly reproduce the software
environment used when these simulations where originally run and may
result in discrepancies due to differences in for instance operating
systems, compilers, and BLAS/LAPACK implementations.

## R Package

If you want to build and experiment with the package, you can do so by
calling

     R CMD INSTALL  .

## Data

The data sets used for the project are not stored on this repository and
have to be downloaded by running the script found in
[`data-raw/`](data-raw/). This does not apply when you use the
singularity container, however, since the data sets are stored inside it
(and could technically be retrieved from it too).
