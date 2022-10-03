# Introduction of Alpha and Delta SARS-CoV-2 VoCs into Switzerland

This is a repository for code used in the manuscript "Importation of Alpha and Delta variants during the SARS-CoV-2 epidemic in Switzerland: phylogenetic analysis and intervention scenarios" by Martina L Reichmuth, Emma B Hodcroft, Julien Riou, and Christian L Althaus.

### How this code is used (see paper for more description):

1. All Alpha (Nextstrain clade 20I) and Delta (Nextstrain clades 21A, 21I, and 21J) sequences sampled in Switzerland prior to 31 Mar 2021 and 31 July 2021, respectively, were selected, using [allClusterDynamics_faster.py](scripts/allClusterDynamics_faster.py) and [getClusterFiles.r](scripts/getClusterFiles.r).
Note that allClusterDynamics.faster.py is from [CoVariants.org](https://covariants.org/) and an updated version can be found [at the CoVaraints github](https://github.com/hodcroftlab/covariants/blob/master/scripts/allClusterDynamics_faster.py).

2. These sequences were then run through the Nextstrain ncov pipeline, [as modified for CoVariants.org](https://github.com/emmahodcroft/ncov_2021), to generate “focal” phylogenies. 
Up to 10,000 global sequences that are most genetically similar to the Swiss sequences were included using the Nextstrain code for [proximity](https://github.com/nextstrain/ncov/blob/master/scripts/get_distance_to_focal_set.py) & [priorities](https://github.com/nextstrain/ncov/blob/master/scripts/priorities.py). In the builds, only samples prior to 31 Mar 2021 or 31 July 2021 were included for Alpha and Delta builds, respectively.

3. The resulting trees are read in and 'collapsed' [as described previously](https://www.nature.com/articles/s41586-021-03677-y#Sec8), in that subtrees that contain only sequences from a single country are collapsed into the parental node recursively. This results in what we refer to internally as 'pies' where a node is respresented by the proportion of sequences that have collapsed into it from each country. This is done by [tree_pie_plot.py](scripts/tree_pie_plot.py).

4. The number of estimated introductions of Alpha & Delta into Switzerland is calculated based on a 'liberal' and 'conservative' approach (see paper for details) in [analyze_slices.py](scripts/analyze_slices.py). The number and dates of introductions are then fed into the transmission model.


For Emma, the original repo (which includes data that can't be publicly shared) is [here](https://github.com/emmahodcroft/2021_Delta). The original build auspice file names are ncov_20I.Alpha.V1-swiss-2021-03-31.json and ncov_21A.Delta-swiss-2021-07-31.json.
