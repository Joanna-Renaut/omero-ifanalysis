# IF OMERO Analysis Pipeline

Create your virtual environment e.g for mac:

    python<version> -m venv venv

    source venv/bin/activate

Install:

    pip install requirements.txt
    
Create 3 directories, one for each biological repeat call these:
    
    repeat-1

    repeat-2

    repeat-3

place your Omero final data files in the appropriate folder.

In the code on:
- line 117, enter your control genotype e.g. 'WT'
- on 118 enter your knockout genotype(s) e.g. ['BRCA', 'BRAF']
- on 119 enter your main control e.g. 'NT' 
- on 120 enter all your controls in the order you would like them to appear in figures. 
- Finally on 121 enter the threshold you'd like your plot bars to change colour to signify a good ratio of knockdown, 
e.g. 0.7

Run the analysis.

### Results output
You will have:

- a ratio figure per KO genotype
- a final csv of the normalised absolute cell count per condition / per genotype with the average ratio across 
- repeats and the sem.
- a p-values file showing the statistics using a Welsch's t-test.
