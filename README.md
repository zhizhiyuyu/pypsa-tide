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

## About

### Built With

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
   conda env create -f environment.yaml
```
3. install edited atlite:
```bash
   cd ./atlite
   pip install -e .
```

4. install package
```bash
  #uptide 
  pip install uptide

  # the Aviso FES package
  conda install -c anaconda cmake
  pip install packaging
  pip install wheel
  pip install --no-build-isolation git+https://github.com/CNES/aviso-fes/

  # solver  
  conda install -c gurobi gurobi
  conda install -c conda-forge ipopt coincbc
  conda install glpk  
```
5. download data 
- [ocean_tide_extrapolated](https://imperiallondon-my.sharepoint.com/:f:/r/personal/zz5322_ic_ac_uk/Documents/ocean_tide_extrapolated?csf=1&web=1&e=jyCYu1)
  
## Instruction

### atlite
```bash
# create cutout
# fes_path is the path to ocean_tide_extrapolated filefold
import atlite
cutout = atlite.Cutout(path="test_tide.nc",
                       module=["tide"],
                       fes_path = 'path_to_ocean_tide_extrapolated',
                       x = slice(-13.6913, 1.7712), 
                       y = slice(49.9096, 60.8479),
                       time="2011-01-01", 
                       dt='h' 
                       )
cutout.prepare()

# calculate maximum install capacity (MWh/km²) and potential energy extracted from tidal range (kWh/m²)
capacity_per_sqkm = cutout.tidalrange(return_capacity=True) 
potential_energy_per_sqkm = cutout.tidalrange(return_capacity=False) 

# sum of generation
generation = cutout.tidalplant(lagoon=’20MW’, operation=’two-way’, area= available_area)

# static capacity factor
cf = cutout.tidalplant(lagoon=’20MW’, operation=’two-way’, area= available_area , capacity_factor=True)

# generation time series and
profiles, capacities = cutout.tidalplant(lagoon=’20MW’, operation=’two-way’, area= available_area , layout=layout , matrix=availability , return_capacity= True)
 
# cost calculation (eur/kW)
cost = atlite.tide.calculate_cost(capacity, area)

```
### pypsa-eur

```bash
conda activate pypsa-eur-edi
snakemake -call solve_elec_networks --configfile path-to-configfile
# dry run
# snakemake -call solve_elec_networks --configfile path-to-configfile -n
```
For the configuration file in config/renewable, please run with config.renewable2019.yaml first


