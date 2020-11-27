# Contig Analysis Scripts
Several small python scripts to analyse long read and short read data for scaffolding MHC contigs

## Installation
The scripts have some dependencies on certain Python libraries. If you use conda you can simply install these dependencies in a new conda environment with the provided list of necessary packages.
```bash
conda create --name contig_analysis --file contig_analysis_env.yml
conda activate contig_analysis
```
The scripts do not require many libraries. The main one being biopython, so manual installation with pip should not be too complicated.
When you have all the dependencies, you should be able to simply clone this repository and run the scripts as detailed below.
```bash
git clone https://github.com/DiltheyLab/ContigAnalysisScripts
cd ContigAnalysisScripts
```

## Data generation
These scripts require [PAF formatted ](https://github.com/lh3/miniasm/blob/master/PAF.md) input. This is the default output format of minimap2. 
The scripts' are intended for the purpose of visualizing a region in a larger genome that has been analyzed beforehand, such that a number of high-confidence contigs and a set of Nanopore long reads exist. In the following the high-confidence contiguous sequences will be referred to simply as "contigs" and the file containing them as `contigs.fasta`. The long reads libraries will be referred to as `lr_lib1.fastq`, `lr_lib2.fastq` etc. .
For mapping your data it is advised to add the contigs to the reference genome and mask the parts of the reference with N's that are very similar to your contigs, i.e. the part of the reference that you want to reassemble.
```bash
minimap2 -x map-ont reference_maskedSimilarRegion_contigsAdded.fasta lr_lib1.fastq > lr_lib1.paf
```

## paf2svg.py
This script can be used to visualize all long reads where specific contigs map to. It is advised to use the `--whitelist` parameter to filter for only a few contigs, otherwise the image can be become very big. The `--alignreads` parameter should be used to perform fast pseudo-alignment which should give reads an offset such that the contigs on these reads are aligned nicely.
Generate a `contig_whitelist.txt`file containing contigs (one per line) that you would like to analyse like the following (e.g. with linename SSTO)
```bash
114SSTO
22SSTO
```
You can speficify several different libraries (in .paf format) to analyze them at the same time. You also have to provide the file with the contigs and a linestring. The linestring is stripped from the right of the contig-name and the character before the linestring is used for the coloring of the contig in the generated image. The following example demonstrates this:

(e.g. with linename SSTO)
- 10SSTO -> "0" is used for coloring
- 34SSTO -> "4" is used for coloring

Run the following command to generate a [SVG](https://de.wikipedia.org/wiki/Scalable_Vector_Graphics) image, that can be viewed and edited with Inkscape or Adobe Illustrator.
```bash
python erate2SVG.py --whitelist contig_whitelist.txt --alignreads lr_lib1.paf lr_lib2.paf lr_lib3.paf contigs.fasta linename output.svg
```


![Example image of erate2SVG](https://github.com/TorHou/ContigAnalysisScripts/blob/master/example_erate2SVG.png)

## less_naive_scaffolding

![Example image of less_naive_scaffolding](https://github.com/TorHou/ContigAnalysisScripts/blob/master/example_less_naive_scaffolding.png)
