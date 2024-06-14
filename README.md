# bovine_lavage_scrna

Analysis code for the bovine lavage project.
Documentation will be added as the project progresses.

## Guidelines
- Branches will be named scripts-initals (i.e. scripts-da | scripts-gg) unless this does not work well
- Work on the analysis then submit a pull request when you are ready to merge with main (don't worry about this yet)
- Do not commit directly to main (exception would be correcting typos etc; don't worry about this yet)

<br>

## File structure:
- [:file\_folder: metaData](/metaData) contains relevant metadata files including cell type annotations and alternate sample names for analysis
- [:file\_folder: main](/) contains the scripts used in the analysis. One file for each major analysis step

<br>

## Setting up Git

```sh
#navigate to working directory
cd /pl/active/dow_lab/dylan/bov_lav_scRNA/

#clone the repo -- only need to do this once (can use git fetch/pull once cloned)
git clone https://github.com/dyammons/bovine_lavage_scrna.git

#rename the repository
mv bovine_lavage_scrna scripts-gg

#enter repo
cd scripts-gg/

#create your own branch to work on!
git checkout -b gg
```

<br>

## Using the singularity contianer
To run an interactive R session:
```sh
singularity run -B $PWD/../../../ --fakeroot ../software/r4.3.2-seurat5
```

To submit a job that runs a R script:
```sh
singularity exec -B $PWD/../../../ --fakeroot ../software/r4.3.2-seurat5 Rscript script.R #change script.R to name of script to run
```

<br>

## Directory structure
The project directory is located on the `pl` at `/pl/active/dow_lab/dylan/bov_lav_scRNA`.
The file structure is as follows:

```sh
bov_lav_scRNA
├── external_data
│   └── any data from published sources
├── input
│   └── cellranger output count matrices
├── output
│   └── all output files from analysis
├── scripts-initials
│   ├── .git files
│   ├── logFiles
│   ├── metaData
│   └── where is repository should be cloned!
└── software
    └── singularity sandbox (and other containers)
```
