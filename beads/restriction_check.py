# Checks whether a restriction site is found in any of the bead barcodes

import sys
from itertools import product

comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
site = sys.argv[1]
revcomp_site = "".join(comp[base] for base in site[::-1])
subbarcodes = ['ACGT', 'ACTG', 'AGCT', 'AGTC', 'ATCG', 'CAGT', 'CGAT', 'CGTA', 'CTGA', 'GACT', 'GATC', 'GCTA', 'GTAC', 'GTCA', 'TACG', 'TAGC', 'TCAG', 'TCGA', 'TGAC', 'TGCA']
barcodes = ["".join(i) for i in product(subbarcodes, repeat=4)]
cut = [site in barcode or revcomp_site in barcode for barcode in barcodes]
print(f"{sum(cut)} barcodes cut")

# Test with ACGTACTG and CAGTACGT. Both should be found in cut[1].
# CAGTACGT only found if the reverse complement search works
# Also try GGGGGGGG as negative control
# print(cut[1])
