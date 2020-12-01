import pandas as pd
from pathlib import Path

def calculate_annual_investment(name, r, fn):
    """Calculate the annual investment for an installation given a selected WACC.
    
    Annual investment for the selected installation 'name' is calculated
    using the given WACC as 'discountrate' and properties for the investment
    read from file (CAPEX, FOM in %, lifetime n in years).
    
    Parameters:
    -----------
    name : str
        of the installation investment, e.g. "onwind".
        Requires corresponding entry in file 'fn'.
    r : float
        discount rate or WACC used to calculate the annuity of the investment.
        In %.
    fn : str or pathlib.Path
        Filename or path to a pandas-readable csv containing the cost and investment
        specific information (investment, FOM, lifetime) for the installation.
    
    """
    
    fn = Path(fn)
    assert fn.exists()
    
    costs = pd.read_csv(fn)

    costs = costs[costs["technology"] == name]
    costs = costs.set_index("parameter")

    r = r/100.

    annuity_factor = r/(1. - 1./(r+1.)**(costs.loc["lifetime", "value"]))

    return (annuity_factor + costs.loc["FOM", "value"]/100.)*costs.loc["investment", "value"]