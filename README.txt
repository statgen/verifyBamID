verifyBamID 1.1.3 -- verify identity and purity of sequence data
(c) 2010-2014 Hyun Min Kang, Goo Jun, and Goncalo Abecasis

Forked by simon.white@helix.com  March 14th 2017

The code has been forked in order to produce an additional output file <prefix.id>.
The extra file is used by the helix r2v pipeline to create a identity check to ensure merged lanes
belong to the same individual.

The format of the tab delimited .id file is as follows:

#CHROM	POS	A1	A2	AF	#REF	#ALT	#OTHERS

Where:
A1       = allele 1 (ref)
A2       = allele 2
AF       = allele freq supplied in vcf (thousand genomes)
#REF     = number of reads with A1 allele
#ALT     = number of reads with A2 allele
#OTHERS  = number of reads with other alleles