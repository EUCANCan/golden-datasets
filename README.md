# golden-datasets

Download, extract and format golden datasets to benchmark pipelines

## Requirements

* GNU make
* [miniconda3](https://docs.conda.io/en/latest/miniconda.html)

## Usage

Every steps described below are launched using the  GNU `make` command line tool. To process a command in parallel, the `-j N` option can be used where `N` is the is the maximum number of jobs that `make` will run in parallel.

```shell
make COMMAND -j NCPUS [VARS]
```

| Command | Description |
| --- | --- |
|`download`| Download datasets in `data/datasets.tsv` |
|`bam2fastq`| Convert downloaded `bam` files into `fastq` files |
|`compress`| Compress downloaded `fastq` files with `pigz` |
|`clean`| Remove downloaded data |

## Tutorial

### Download

To download datasets listed in `data/datasets.tsv`, use the `make` command described below. Here we set to run the download process in parallel 4 at a time.

```shell
make download -j 4 INSTALL_DIR=$HOME/golden-datasets
```

By default, the installation directory is set to the `data` folder in the main repository aside the dataset `tsv` files. This behaviour can be change like in the example above where the installation directory is set from the environment variable `INSTALL_DIR`. Environment variables can be passed to the `make`command line tool as positional arguments. 

<aside class="notice">
Several files in those datasets require a synapse account in order to download them. Go to https://www.synapse.org/ in order to create an account. Those credentials will be asked the first time the make command is launched.</aside>


### Convert and compress downloaded files

Mapping files located in the `INSTALL_DIR` folder can be converted to fastq with the command below. The tool `Picard SamToFastq` is used for the conversion.

```shell
make bam2fastq -j 4
```

After the download and conversion, `FASTQ` files can be compressed with `pigz` using the command described below

```shell
make compress -j 4
```

## License

Licensed under the
[GNU General Public License v3.0](https://github.com/EUCANCan/golden-datasets/blob/master/LICENSE) License.
