### get_competing_splice_sites.py  

Identifies splice sites that are both unique to specific pairwise alternative splicing event isoforms and in competition with another splice site for a particular partner.  Extracts and scores splice site sequences and dinucleotides, and compares the scores of competitors.  If input data from more than one species is provided, information on the splice sites from all species will be reported, along with a comparison of the scores.


### prep_splice_site_metaplots.py

Wraps bedtools to use a regions bed file (consisting of regions with the same length) to extract data from a bedGraph file (e.g. phyloP scores) in a form convenient for plotting.  
Note that this will likely be replaced with an alternative that wraps bwtool (by Andy Pohl - see https://github.com/CRG-Barcelona/bwtool/ and http://bioinformatics.oxfordjournals.org/content/early/2014/02/16/bioinformatics.btu056) which seems much faster for this application and generally excellent.
