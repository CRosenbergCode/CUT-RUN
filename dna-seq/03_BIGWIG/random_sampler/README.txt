This is a simple tool I made for generating random bed files for a genome. This was used for debugging issues we were having with our bigwigs, generating DeepTools plots.

Essentially, all it does is create a bed file with random regions of a given size, so that you can use that bed file when generating DeepTools heatmaps/profiles and see how two samples compare to one another in overall RPGC, CPM, etc. This is purely for debugging.

Run:

python pick_regions.py -h

to see how to use it. Example region files are available in this folder. Here is an example of using it, using the example region files. This works for Aedes aegypti.

python pick_regions.py --contigsizes contig_lengths.txt --contignames contig_names.txt --range 4000 -n 10000 -o randomSamples.bed

This will create a file, randomSamples.bed, with 10000 total randomly-picked regions. Each region is 4000 base pairs long. The contig_lengths.txt file contains the lengths of each contig, and the contig_names.txt file contains the names of the corresponding contigs in the same order. This script could later be updated to accept a single combined file with this data.