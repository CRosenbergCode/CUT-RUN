To create Aedes, the VectorBase-68_AaegyptiLVP_AGWG.gff file (downloaded from Vectorbase and available in the NAS) was converted to a gtf file using AGAT. Use of AGAT is important as other tools generally do not preserve UTRs.

Ex: agat_convert_sp_gff2gtf.pl --gff VectorBase-68_AaegyptiLVP_AGWG.gff  -o Aedes68.gtf

To create the AedesGenes.gtf file, grep was used to select only lines with genes from the Aedes68.gtf file.

Ex: grep protein_coding_gene Aedes68.gtf > AedesGeneExtras.gtf

To create AedesGenesExtras.gtf, the method above was used and also applied to long non-coding RNAs (lncs) and pseudogenes.

Ex:
grep protein_coding_gene Aedes68.gtf > AedesGeneExtras.gtf
grep lncRNA Aedes68.gtf >> Aedes68GenesExtra.gtf
grep "\spseudogene\s" Aedes68.gtf >> Aedes68GenesExtra.gtf

To create Aedes68Exon3Prime.gtf, the method above was used to extract exons and 3’ UTRs

Ex:

grep 'three_prime_UTR' Aedes68.gtf >> Aedes68Exon3Prime.gtf
grep 'three_prime_UTR' Aedes68.gtf >> Aedes68Exon3Prime.gtf
