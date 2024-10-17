#%%
from pathlib import Path

import pandas as pd

from loopy.drivers.spaceranger import run_spaceranger
from loopy.feature import FeatureAndGroup

#%%
directory = Path("/Users/jianingyao/Desktop/Research/Biostatistics_JHU/PhD/Data/samui-data/")
out = Path("/Users/jianingyao/Desktop/Research/Biostatistics_JHU/PhD/Data/samui-data/samui-output")
samples = ["Visium_CKp25_rep1", "Visium_CKp25_rep2", "Visium_CKp25_rep3", "Visium_CKp25_rep4"]
samples = ["Visium_CKp25_rep4"]

channels = [
    "DAPI",
    "channel2",
    "Protein",
]

out = Path("/Users/jianingyao/Desktop/Research/Biostatistics_JHU/PhD/Data/samui-data/samui-output")

# importantFeatures = (
#     [
#         FeatureAndGroup(feature="GFAP", group="genes"),
#         FeatureAndGroup(feature="OLIG2", group="genes"),
#         FeatureAndGroup(feature="TMEM119", group="genes"),
#         FeatureAndGroup(feature="RBFOX3", group="genes"),
#     ],
# )

for name in samples:
    s = (
        run_spaceranger(
            name,
            directory / name,
            out,
            tif=directory / (name + ".tif"),
            channels=channels,
            defaultChannels=dict(blue="DAPI", green="GFAP", red="NeuN"),
        )
        .add_coords(pd.read_csv(directory / name / "cellCoords.csv"), name="cells", mPerPx=0.497e-6)
        .add_csv_feature(
            pd.read_csv(directory / name / f"{name}_filtered.csv"),
            name="CellsFiltered",
            coordName="cells",
            dataType="categorical",
        )
        .add_csv_feature(
            pd.read_csv(directory / name / f"{name}_unfiltered.csv"),
            name="CellsUnfiltered",
            coordName="cells",
            dataType="categorical",
        )
        .set_default_feature(group="genes", feature="GFAP")
        .write()
    )