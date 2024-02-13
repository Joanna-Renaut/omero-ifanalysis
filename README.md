# IF OMERO Analysis Pipeline

### IFAnalysis

First get the IFAnalysis pipeline by Helfrid Hocchegger: https://github.com/HocheggerLab/IFanalysis
Install via

    pip install -e <path to repo on computer>

Create 3 directories, one for each biological repeat call these:
    
`repeat-1`

`repeat-2`

`repeat-3`

place your Omero final data files in the appropriate folder.

Run the analysis.

You will have multiple figures, a final csv of the abs cell counts and relative cell counts (ratios) in the final-data 
folder. A final stats file with the p-values for each repeat and a csv for means and std deviations for abs 
cell count and relative cell counts (ratios) in the stats folder.