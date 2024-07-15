
import functools
import logging
import time

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import configure_logging
from dask.distributed import Client
from pypsa.geo import haversine
from shapely.geometry import LineString

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_tidal_profiles", technology="tidal")
    configure_logging(snakemake)

    nprocesses = int(snakemake.threads)
    noprogress = snakemake.config["run"].get("disable_progressbar", True)
    noprogress = noprogress or not snakemake.config["atlite"]["show_progress"]
    params = snakemake.params.tidal
    correction_factor = params.get("correction_factor", 1.0)
    scaling_factor =  params.get("scaling_factor", 1.0)
    resource = params["resource"] 
    if correction_factor != 1.0:
        logger.info(f"correction_factor is set as {correction_factor}")

    if nprocesses > 1:
        client = Client(n_workers=nprocesses, threads_per_worker=1)
    else:
        client = None

    sns = pd.date_range(freq="h", **snakemake.config["snapshots"])
    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=sns)
    regions = gpd.read_file(snakemake.input.regions)
    assert not regions.empty, (
        f"List of regions in {snakemake.input.regions} is empty, please "
        "disable the corresponding renewable technology"
    )
    # do not pull up, set_index does not work if geo dataframe is empty
    regions = regions.set_index("name").rename_axis("bus")
    buses = regions.index
    excluder = atlite.ExclusionContainer(crs=3035, res=100)
   
    kwargs = dict(nprocesses=nprocesses, disable_progressbar=noprogress)
    if noprogress:
        logger.info("Calculate landuse availabilities...")
        start = time.time()
        availability = cutout.availabilitymatrix(regions, excluder, **kwargs)
        duration = time.time() - start
        logger.info(f"Completed availability calculation ({duration:2.2f}s)")
    else:
        availability = cutout.availabilitymatrix(regions, excluder, **kwargs)

    func = getattr(cutout, resource.pop("method"))
    if client is not None:
        resource["dask_kwargs"] = {"scheduler": client}
    capacity_factor = correction_factor * func(capacity_factor=True,**resource)

    data = pd.read_excel(resource["path"]).rename(columns={'Longitude':'x', 'Latitude':'y'})
    layout = cutout.layout_from_capacity_list(data, col = "Project Capacity")

    profile, capacities = func(
        matrix=availability.stack(spatial=["y", "x"]),
        layout = capacity_factor * layout,
        index=buses,
        per_unit=True,
        return_capacity=True, 
        **resource,
    )

    p_nom_max = availability @ layout

    logger.info("Calculate average distances.")
    layoutmatrix = (layout * availability).stack(spatial=["y", "x"])
    
    coords = cutout.grid[['x', 'y']]
    bus_coords = regions[['x', 'y']]

    average_distance = []
    centre_of_mass = []
    for bus in buses:
        row = layoutmatrix.sel(bus=bus).data
        nz_b = row != 0
        row = row[nz_b]
        co = coords[nz_b]
        distances = haversine(bus_coords.loc[bus],  co)
        average_distance.append((distances * (row / row.sum())).sum())
        centre_of_mass.append(co.values.T @ (row / row.sum()))

    average_distance = xr.DataArray(average_distance, [buses])
    centre_of_mass = xr.DataArray(centre_of_mass, [buses, ('spatial', ['x', 'y'])])


    ds = xr.merge([(correction_factor * profile).rename('profile'),
                    capacities.rename('weight'),
                    p_nom_max.rename('p_nom_max'),
                    capacity_factor.rename('capacity_factor'),
                    average_distance.rename('average_distance')])



    logger.info("Calculate underwater fraction of connections.")
    offshore_shape = gpd.read_file(snakemake.input["offshore_shapes"]).unary_union
    underwater_fraction = []
    for bus in buses:
        p = centre_of_mass.sel(bus=bus).data
        line = LineString([p, regions.loc[bus, ["x", "y"]]])
        frac = line.intersection(offshore_shape).length / line.length
        underwater_fraction.append(frac)

    ds["underwater_fraction"] = xr.DataArray(underwater_fraction, [buses])

    # select only buses with some capacity and minimal capacity factor
    ds = ds.sel(
        bus=(
            (ds["profile"].mean("time") > params.get("min_p_max_pu", 0.0))
            & (ds["p_nom_max"] > params.get("min_p_nom_max", 0.0))
        )
    )
    if "clip_p_max_pu" in params:
        min_p_max_pu = params["clip_p_max_pu"]
        ds["profile"] = ds["profile"].where(ds["profile"] >= min_p_max_pu, 0)

    ds.to_netcdf(snakemake.output.profile)

    if client is not None:
        client.shutdown()
