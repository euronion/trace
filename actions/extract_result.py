# SPDX-FileCopyrightText: 2023 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

import re

import numpy as np
import pandas as pd
import pypsa

scenario = snakemake.params["scenario"]

l = []

n = pypsa.Network()
n.import_from_netcdf(snakemake.input["network"])

res = n.generators_t["p"].filter(regex="wind|pvplant")
res = res.loc[:, (res != 0.0).any()]
res_opt = res.divide(n.generators["p_nom_opt"])
res_opt = res_opt.loc[:, (res_opt != 0.0).any(axis=0)].dropna(axis=1)


# In[ ]:


# Use unit of load bus and add "h" or remove "/h".
# Bus units should be in flow units
units = n.buses.loc[n.loads.bus.unique()]["unit"].unique()
len(units) == 1, "Only exactly one unit for loads currently supported"
load_unit = units[0]
if load_unit.endswith("/h"):
    load_unit = load_unit.replace("/h", "")
else:
    load_unit += "h"


# In[ ]:


## General statistics
category = "general"

l.append([category, "Total system cost", n.objective])

# Price per unit of load delivered
l.append(
    [
        category,
        f"Cost per {load_unit} delivered",
        n.objective / n.loads_t["p"].sum().sum(),
    ]
)

# Total demand
l.append(
    [
        category,
        "Total demand",
        n.loads_t["p"].sum().sum(),
    ]
)


# In[ ]:


## RES
category = "RES"

# Installed RES capacity
l.append(
    [
        category,
        "Total installed capacity",
        n.generators["p_nom_opt"].sum(),
    ]
)

# Installed RES capacities (primary energy: electricity or heat from csp)
for g, v in n.generators.groupby("carrier")["p_nom_opt"].sum().items():
    l.append([category, f"Total installed capacity {g}", v])

# Generated electricity from CSP power block or conventional renewables
(
    n.links_t["p1"]
    .filter(like="csp-tower power block", axis="columns")
    .abs()
    .sum()
    .sum()
    + n.generators_t["p"].filter(regex="solar-utility|wind", axis="columns").sum().sum()
)

# Generated electricity from CSP power block or conventional renewables
elec = pd.concat(
    [
        n.links_t["p1"]
        .filter(like="csp-tower power block", axis="columns")
        .abs()
        .sum(),
        n.generators_t["p"].filter(regex="solar-utility|wind", axis="columns").sum(),
    ]
)

elec = elec.groupby(lambda x: re.sub("\s\d+\s?", "", x)).sum()

# Produced electricity
l.append(
    [
        category,
        "Total generated electricity",
        elec.sum(),
    ]
)

# Energy surplus factor
l.append(
    [
        "general",
        "Energy surplus factor",
        elec.sum() / n.loads_t["p"].sum().sum(),
    ]
)

for idx, v in elec.items():
    # Generated electricity from all sources
    l.append([category, f"Total generated electricity from {idx}", v])

# Curtailed RES energy
l.append(
    [
        category,
        "Curtailed RES energy",
        (n.generators_t["p_max_pu"] * n.generators["p_nom_opt"] - n.generators_t["p"])
        .sum(axis=1)
        .sum(),
    ]
)


# In[ ]:


# costs
category = "cost"

t = getattr(n, "generators")
t_t = getattr(n, "generators_t")
t = (
    (
        (t["p_nom_opt"].clip(lower=0.0) * t["capital_cost"])
        + (t_t["p"] * t["marginal_cost"]).sum()
    )
    .replace(0.0, np.nan)
    .dropna()
)

t = list(zip([category] * len(t), t.index, t.to_list()))
l.extend(t)

t = getattr(n, "links")
t_t = getattr(n, "links_t")
t = (
    (
        t["p_nom_opt"].clip(lower=0) * t["capital_cost"]
        + (t_t["p0"] * t["marginal_cost"]).sum()
    )
    .replace(0.0, np.nan)
    .dropna()
)

t = list(zip([category] * len(t), t.index, t.to_list()))

l.extend(t)

t = getattr(n, "stores")
t_t = getattr(n, "stores_t")
t = (
    (
        t["e_nom_opt"].clip(lower=0) * t["capital_cost"]
        + (t_t["e"] * t["marginal_cost"]).sum()
    )
    .replace(0.0, np.nan)
    .dropna()
)
t = t.groupby(
    lambda x: re.sub(
        "\sconvoy \d+ .*", "", x.replace(" (exp)", "").replace(" (imp)", "")
    )
).sum()

t = list(zip([category] * len(t), t.index, t.to_list()))

l.extend(t)


# In[ ]:


# capacity factors

# Annual CF
t = res.to_numpy().sum() / (
    n.generators.filter(regex="wind|pvplant|csp-tower", axis=0)["p_nom_opt"].sum()
    * n.snapshots.shape[0]
)

l.append(["capacity factor", "RES annual CF", t])

# Average annual CF
t = res_opt.to_numpy().mean()

l.append(["capacity factor", "RES average CF", t])

c = "stores"
t = getattr(n, c)
t_t = getattr(n, c + "_t")

# Cleanup rounding errors and small values != 0
t = t["e_nom_opt"].round(decimals=6)
t_t = t_t["e"].round(decimals=6)

t = t_t.divide(t).dropna(axis=1).mean()
t = list(zip(["capacity factor"] * len(t), t.index, t.to_list()))

l.extend(t)

c = "links"
t = getattr(n, c)
t_t = getattr(n, c + "_t")

# Cleanup rounding errors and small values != 0
t = t["p_nom_opt"].round(decimals=6)
t_t = t_t["p0"].round(decimals=6)

t = t_t.divide(t).dropna(axis=1).mean()
t = list(zip(["capacity factor"] * len(t), t.index, t.to_list()))

l.extend(t)

c = "generators"
t = getattr(n, c)
t_t = getattr(n, c + "_t")

t = t_t["p"].divide(t["p_nom_opt"].clip(lower=0.0)).dropna(axis=1).mean()
t = list(zip(["capacity factor"] * len(t), t.index, t.to_list()))

l.extend(t)


# In[ ]:


# build capacities


# In[ ]:


c = "generators"
t = getattr(n, c)

t = t["p_nom_opt"].clip(lower=0.0).replace(0.0, np.nan).dropna()
t = list(zip(["installed capacity"] * len(t), t.index, t.to_list()))

l.extend(t)

c = "links"
t = getattr(n, c)
t = t["p_nom_opt"].clip(lower=0.0).replace(0.0, np.nan).dropna()

t = list(zip(["installed capacity"] * len(t), t.index, t.to_list()))

l.extend(t)

c = "stores"
t = getattr(n, c)
t = t["e_nom_opt"].clip(lower=0.0).replace(0.0, np.nan).dropna()

t = list(zip(["installed capacity"] * len(t), t.index, t.to_list()))

l.extend(t)


# In[ ]:


# List to dataframe

df = pd.DataFrame(l, columns=["category", "subcategory", "value"])


# In[ ]:


# Add experiment specific information
for k, v in snakemake.wildcards.items():
    df[k] = v
df["wacc"] = scenario["wacc"]

df = df.set_index(
    [
        "scenario",
        "year",
        "esc",
        "exporter",
        "importer",
        "wacc",
        "category",
        "subcategory",
    ],
    verify_integrity=True,
)


# In[ ]:


# Dataframe to file

df.to_csv(snakemake.output["results"], sep=";")


# In[ ]:


import seaborn as sns

sns.set_style("whitegrid")
import matplotlib.pyplot as plt

# In[ ]:


# min electricity generation
t = res.sum(axis=1).min()

l.append(["general", "minimum RES generation", t])

plt.figure(figsize=(30, 30))
ax = plt.gca()

res = n.generators_t["p"].filter(regex="pvplant|wind", axis=1)
res = res.loc[:, (res != 0.0).any()]

real = res.sum(axis=1)
theo_max = (
    n.generators_t["p_max_pu"]
    .multiply(n.generators["p_nom_opt"])[res.columns]
    .sum(axis=1)
)
diff = real - theo_max

theo_max.name = "Available generation"
real.name = "Realised generation"
diff.name = "Curtailed"

theo_max.plot(ax=ax)
real.plot(ax=ax)
diff.plot(ax=ax)

plt.legend()


# In[ ]:


if "shipping" in n.name:
    n.links_t["p0"].filter(regex="ship convoy [0-9]* loading").sum(axis=1).plot(
        figsize=(30, 10), title="ship loading"
    )


# In[ ]:


if "shipping" in n.name:
    n.links_t["p0"].filter(regex="ship convoy [0-9]* unloading").sum(axis=1).plot(
        figsize=(30, 10), title="ship unloading"
    )


# In[ ]:


if "shipping" in n.name:
    n.stores_t["e"].filter(like="cargo").sum(axis=1).plot(
        figsize=(30, 10), title="ship cargo"
    )


# In[ ]:


cs = {c for c in n.links.index if "ship convoy" not in c}


# In[ ]:


c1 = {c for c in cs if c.startswith("battery inverter")}
c2 = {c for c in cs if c.startswith("seawater desalination")}
c3 = cs - c1 - c2

c1 = list(c1)
c2 = list(c2)
c3 = list(c3)


# In[ ]:


n.links_t["p0"].index = pd.to_datetime(n.links_t["p0"].index)


# In[ ]:


for c in [c1, c2, c3]:
    if c:
        n.links_t["p0"][c].plot(figsize=(50, 20))
        plt.ylim(
            0,
        )


# In[ ]:


cs = {c for c in n.stores_t["e"].columns if "cargo" not in c}


# In[ ]:


c1 = {
    c
    for c in cs
    if c.startswith("hydrogen storage tank")
    or c.startswith("Buffer: hydrogen (g) demand")
}
c2 = {
    c
    for c in cs
    if c.startswith("clean water tank storage") or c.startswith("battery storage")
}
c3 = cs - c1 - c2

c1 = list(c1)
c2 = list(c2)
c3 = list(c3)


# In[ ]:


for c in [c1, c2, c3]:
    if c:
        n.stores_t["e"][c].plot(figsize=(50, 20), title="Stores")
        plt.ylim(
            0,
        )
