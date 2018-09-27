## Edchecker
The name does not matter right now :)

### Summary:

TL,DR: See [dataset/editing_map](dataset/editing_map) for the dataset !

### Content

This repository only contains code for creating a dataset of known RNA editing sites in plant mitochondrial genomes, and a few utility function that I need to encode a MSA as a matrix.
In the next few hours (when I am free), I will start adding implementation of the methods described in the document. 
<pre>
.
+-- setup.py     ==> not really useful anymore as I am using click 
+-- README.md    ==> this file
+-- .gitignore 
+-- misc
|   +-- genelist  ==> list of accepted gene list 
|   +-- dbraw/editing.txt     ==> content raw database file, 
|                                 essentially a text dumb of the dead REDidb database
+-- dataset
|   +-- Edit.db      ==> sqlite database that contains a comprehensible list of 
|                        RNA editing in various genomes (see details below)
|   +-- prot_align/  ==> folder that contains all proteins MSA (obtained with mafft) 
|                        that I used as template for the codon-by-codon alignment
|   +-- genomes/     ==> Folder that contains raw mitochondrial genomes in genbank format
|
|   +-- <b>editing_map  ==> dataset of interest that give the edited positions in the alignment</b>
+-- AlignFormat.py  ==> Format alignment into matrix
+-- DBManager.py    ==> Parse the text database of RNA editing and save it into an sqlite database
+-- DataParser.py   ==> Given a set of genbank file, parse the file using only 
                        genes of interest and extract all potentially edited site, then save them
</pre>

### Potential issues with the true positive dataset

Despite my best efforts, it is still possible that the true positive dataset might be incomplete or could include a few false positives. I have added safeguards to ensure that only edited positions, where `C` is the nucleotide found in the genomic sequence, are considered.
This issue is caused by the fact that several `misc_feature` (used for additional notes) for RNA editing in the genbank file do not make any senses and most of the time RNA editing is present but not even indicated in the annotation. On the other hand, the RNA editing database I am pulling additional information from is old and use genome version that are deprecated on NCBI.

### Building a new dataset

You can run the script (```DataParser.py```) if you want to build a new dataset (ex, by adding new genomes). All scripts come with a help option.
run ```python SCRIPT_NAME --help```

## C2USimul

C2USimul simulates a dataset (codon alignment) with RNAediting. Codon sequences are simulated using pyvolve and mutation rates provided or drawn from an exisiting alignment. 

An example of command to run C2USimul to return 2 'simulated genes' with RNA editing at ~15% or C positions. Note that for each position of the protein alignment a random clade is selected and subjected to editing. The lower the depth of their corresponding LCA, the higher is the chance of selection. This is to prevent editing at orphan genomes, since it is not what we observe in real dataset (because editing site are often conserved).
```
python C2USimul.py  --gtree misc/tree_simul_example.nw  --gsize 2 --glen_range 100 400 -wd test --delrate 0.3 --protlike --from_al dataset/editing_map/atp1 --eprob 0.15
```
