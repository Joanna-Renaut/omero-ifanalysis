# IF OMERO Analysis Pipeline

### IFAnalysis

First get the IFAnalysis pipeline by Helfrid Hocchegger: https://github.com/HocheggerLab/IFanalysis

Create your virtual environment e.g for mac:

    python<version> -m venv venv

    source venv/bin/activate

Install via

    pip install -e <path to IFAnalysis folder on computer>

Next:

    pip install requirements.txt
    
Create 3 directories, one for each biological repeat call these:
    
    repeat-1

    repeat-2

    repeat-3

place your Omero final data files in the appropriate folder.

In the code on line 65, enter your control conditions in the order you would like them to appear in figures.

Run the analysis.

You will have multiple figures, a final csv of the abs cell counts and relative cell counts (ratios) in the final-data 
folder. A final stats file with the p-values for each repeat and a csv for means and std deviations for abs 
cell count and relative cell counts (ratios) in the stats folder.
