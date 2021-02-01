#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pypsa
import numpy as np
import pickle
from pathlib import Path

from _helpers import configure_logging
import logging
logger = logging.getLogger(__name__)

from pypsa.linopt import define_constraints, linexpr, get_var

def extra_functionality(n, snapshots):
    """... for solving the network; from pypsa-eur-sec."""
    add_battery_constraints(n)
    
    if n.generators.index.str.contains('convoy').any():
        add_shipping_constraints(n)

def add_battery_constraints(n):
    """Constraint for battery inverter capacity.
    
    Assumed battery inverter is bidirectional.
    Enforce same capacity for charging/discharging link.
    
    From pypsa-eur-sec. (2020-12-14)
    """
    
    chargers = network.links.query("name.str.startswith('battery inverter') "
                    "and bus1.str.startswith('battery') "
                    "and p_nom_extendable", engine='python').index

    dischargers = network.links.query("name.str.startswith('battery inverter') "
                        "and bus0.str.startswith('battery') "
                        "and p_nom_extendable", engine='python').index

    link_p_nom = get_var(n, "Link", "p_nom")

    lhs = linexpr((1,link_p_nom[chargers]), (-1, link_p_nom[dischargers].values))
                 
    pypsa.linopt.define_constraints(n, lhs, "=", 0, 'Link', 'charger_ratio')

def add_shipping_constraints(n):
    """Constraint for ensuring shipping to work properly.
    
    Ensures same amounts of energy are transferred via the link pairs emulating
    a shipping convoy along a route.
    Shipping efficiencies is considered while loading the convoys (link from commodity store to berth);
    this makes the constraint here simpler.
    """

    unloading = n.links.filter(regex='transport ship convoy [0-9]* unloading', axis=0).index
    loading = unloading.str.replace('unloading','loading')

    link_p_nom = get_var(n, "Link", 'p_nom')

    lhs = linexpr((1, link_p_nom[loading]), (-1., link_p_nom[unloading].values))

    pypsa.linopt.define_constraints(n, lhs, '=', 0, 'Link', 'shipping_constraints')


# In[ ]:


if __name__ == "__main__":

    configure_logging(snakemake)
    
    solver_options = snakemake.config['solver']
    solver_name = solver_options.pop('name')

    # Load additional components
    with open(snakemake.input["additional_components"], "rb") as f:
        override_component_attrs = pickle.load(f)
    network = pypsa.Network(override_component_attrs=override_component_attrs)

    network.import_from_netcdf(snakemake.input["network"])
    
    network.consistency_check()

    logger.info(f'Solving network using solver options: {solver_options}.')
    logger.info("Starting LOPF.")
    status, termination_condition = network.lopf(snapshots=network.snapshots,
                                                 extra_functionality=extra_functionality,
                                                 pyomo=False,
                                                 solver_name=solver_name,
                                                 solver_options=solver_options,
                                                 solver_logfile=snakemake.log['python'])
    logger.info("End of LOPF.")

    if status == 'ok':
        network.export_to_netcdf(snakemake.output["network"])
    else:
        logger.warning(f'Optimisation ended with unexpected status: {status}. '
                       f'Unoptimised network not saved.')

