# PyPSA-Eur + tidal range energy

## Table of Contents

- [PyPSA-Eur + tidal range energy](#pypsa-eur--tidal-range-energy)
  - [Table of Contents](#table-of-contents)
  - [About](#about)
    - [Built With](#built-with)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
  - [Instruction](#instruction)
    - [atlite](#atlite)
    - [pypsa-eur](#pypsa-eur)
    - [snakemake](#snakemake)

## About

A brief introduction to your project. Explain what it does, why it's useful, and who it's for. Mention any related projects or similar solutions.


### Built With

List the technologies, frameworks, and tools used in your project. For example:

- [uptide](https://github.com/stephankramer/uptide)
- [atlite](https://github.com/PyPSA/atlite)
- [PyPSA-eur](https://github.com/PyPSA/pypsa-eur)


## Getting Started

Provide instructions on how to set up and get your project running locally. Include prerequisites, installation steps, and usage examples.

### Prerequisites

List any software, libraries, or dependencies that users need to have installed before they can use your project.

- [gurobi optimization](https://www.gurobi.com/academia/academic-program-and-licenses/)
- [CDS API](https://cds.climate.copernicus.eu/api-how-to)


### Installation

Detailed installation instructions, including any configuration steps. For example:

1. Clone the repository:

```bash
   git clone https://github.com/zhizhiyuyu/pypsa-tide
```

2. install environment:
```bash
   conda env create -f environment.edi.yaml
```

3. install atlite+tide:
```bash
   cd ./atlite
   pip install -e .
```

## Instruction

### atlite

### pypsa-eur

```bash
conda activate pypsa-eur-edi
snakemake -call solve_elec_networks --configfile path-to-configfile
# dry run
# snakemake -call solve_elec_networks --configfile path-to-configfile -n

```
For the configuration file in config/renewable, please run with config.renewable2019.yaml first

### snakemake

