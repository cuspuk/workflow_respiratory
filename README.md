# Snakemake workflow: Respiratory

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/xsitarcik/respiratory/workflows/Tests/badge.svg?branch=main)](https://github.com/xsitarcik/respiratory/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for analysis of respiratory samples.

## Installing and running

To install the workflow, simply git clone the repository into the path you want:

```bash
git clone git@github.com:xsitarcik/respiratory.git
```

Install the following conda environment:

```bash
mamba create -c conda-forge -c bioconda --name snakemake_respiratory snakemake=7.25 peppy snakemake-wrapper-utils
```

**IMPORTANT**: change the directory to the cloned repository - workflow directory. Every relative path mentioned is relative to this directory.

### Preparing data and configuration

First, prepare your data configuration using PEP file, see [samples.csv](config/pep/samples.csv) as an example. Create new file or update the existing one.

Second, you must configure the `config/config.yaml` file:

- Do not forget to update the `pepfile` path, if you created your own PEP file.
- Provide the path to the viral reference panel, i.e. `reference_panel_dirpath`.
- Update other values as desired.

Third, in the path specified by `kraken_dir` there will be created a new directory named as `kraken_db` where the workflow will download a kraken database. To save both, disk space and time required for downloading, link the existing kraken database here, if possible:

```sh
ln -s {existing_kraken_db} {kraken_dir}/{kraken_tag}
```

### Running the workflow

There are many arguments to use when running a snakemake workflow, see [the documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html). Recommended arguments to use, is `--use-conda` to use conda, and the `--conda-prefix` which is the directory where snakemake will install workflow conda environments. Then for example, the command is:

First, it is advised to dry-run the snakemake workflow using `--dry-run`, i.e.:

```shell
snakemake --cores {THREADS} --use-conda --rerun-incomplete --printshellcmds --dry-run
```

A basic summary of outputs and rules is outputted for you to verify. Then, run snakemake without `--dry-run`

```shell
snakemake --cores {THREADS} --use-conda --rerun-incomplete --printshellcmds
```

### Debugging

- run with `--notemp` to ignore `temp()` in snakemake rules.
- run with `--show-failed-logs` to automatically show failed log messages.

### Issues

- use issues in github to report any problems. This way solving any issue will be stored here for other users or future ourselves.
- Watch the repo to subscribe to pull requests, issues and other repo-related stuff. This way you will be notified of any changes.

### Running using SLURM

It is advised to prepare a different directory for cluster logs, for example we can use `./sbatch_logs`. You must create the directory manually as slurm does not create it automatically:

```sh
mkdir -p sbatch_logs
```

Also you will probably need to set executive permissions to the sbatch monitoring script:

```sh
chmod +x sbatch_status.py
```

Then, the recommended command is:

```shell
snakemake --jobs 40 --rerun-incomplete --use-conda --printshellcmds --cluster "sbatch --partition=gen-compute --cpus-per-task=8 --parsable --output=`pwd`/sbatch_logs/%j.log --error=`pwd`/sbatch_logs/%j.err" --cluster-status `pwd`/sbatch_status.py --keep-going --cluster-cancel scancel --retries 3 --dry-run
```

Do not forget to adjust the numbers accordingly for your cluster configuration, and also to run it without the `--dry-run` option.

Some explanations for the above command:

- `%j` is used to create a log for each cluster job named the same as the job itself.
- `--cpus-per-task=8` it is important to also specify the number of threads in the config.
- Snakemake process must be kept alive during the whole analysis as a master process. Sometimes a SLURM job can fail before sending the failure message to the master process and the job will hang forever. To prevent this, snakemake recommends the usage of [cluster-status argument](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-cluster-status), requiring to pass parsable argument to sbatch.
- Snakemake jobs can be expected to fail such as download jobs, but also they can fail due to memory constraints, i.e. when a job requests more memory than expected due to an unusually large .fastq.gz file, you should pass `--retries {XXX}` option to tell snakemake to retry the job. This is further used for example in picard deduplication where more memory is requested for the subsequent retries (up to the maximum memory specified in the config).

### Results

After running, you can create a [report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html), either `.zip` or `.html`:

```shell
snakemake --report my_first_report.zip
```

## Development

Install `snakemake` with `pre_commit`, for example in the environment `snakemake_dev`:

```shell
mamba create -c conda-forge -c bioconda --name snakemake_dev snakemake=7.25 snakemake-wrapper-utils pre_commit peppy
```

Then set up `pre-commit` in the cloned repository:

```bash
pre-commit install
```

Now before any commit, a defined set of actions will be performed, such as linting, formatting, etc.

### Testing and linting

There is set up an automatic workflow running in the `.tests/` directory. Here you should define a config file and provide any test inputs.
Formatting and linting is also automated for both Python scripts and snakemake rules.

### Conventional PRs

When committing, you must follow the [Conventional Commits spec](https://www.conventionalcommits.org/en/v1.0.0/). Each PR is automatically validated by the GH action.

Further, any push (i.e. after merged PR) to the `main` branch creates in a new PR:

- a new release following the [Semantic Versioning](https://semver.org/)
- an automatic changelog as parsed from the commit history
