# Bachelor thesis

This repo contains the bachor thesis of me!  
It's quite big since I'm lazy and just included all the data in here...  

The repo that contains the written `spam-align` tool is added as a submodule under `alignment-evaluation/spam-align` and is also made available on gitlab [here](https://gitlab.gwdg.de/robinwilliam.hundt/spam-align).  


## Contents
### alignment-evaluation
Contains source code that is used to run `spam-align`,`mafft` and `dialign2-2` on the Balibase 3 dataset and evaluate, aggregate, visualize the results

### datasets
Contains the Balibase3 and IRMBase2 datasets. Currently only Balibase is used but IRMBase might be used in the future.

### evaluation-data
The results of running the different alignment tools on the balibase data and scoring the resulting alignments. This folder is created by the tool in `alignment-evaluation`.

### pattern_sets
The pattern sets that are used for the evaluation of `spam-align`.