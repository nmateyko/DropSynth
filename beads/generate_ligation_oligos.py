# Make oligos for testing split-and-pool ligation synthesis of barcoded beads

import datetime

def revcomp(seq):
  comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  return "".join(comp[base] for base in seq[::-1])

subbarcodes = ['ATCG', 'CGAT', 'GATC', 'AGTC']
anchor = 'TTGTTAGGTCTAGA'
free = 'TGAGACCAACCTC'
end = 'GCTCTTCGCAATC'
date = datetime.date.today()
rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

with open(f'ligation_oligos_{str(date)}.csv', 'w') as f:
  # write header
  f.write('Well Position,Name,Sequence\n')
  # write complementary anchor oligos
  for i, barcode in enumerate(subbarcodes):
      f.write(f'{rows[0]}{i + 1},anchor_rev_{barcode},{revcomp(barcode)}{revcomp(anchor)}\n')
  # write free oligo top strands
  i = 0
  for sb1 in subbarcodes:
    for sb2 in subbarcodes:
      row = rows[i // 8 + 2]
      f.write(f'{row}{i % 8 + 1},free_{sb1}_{sb2},{sb1}{sb2}{free}\n')
      i += 1
  # write free oligo bottom strands
  for i, barcode in enumerate(subbarcodes):
    f.write(f'{rows[5]}{i + 1},free_rev_{barcode},{revcomp(free)}{revcomp(barcode)}\n')
  # write end oligos
  for i, barcode in enumerate(subbarcodes):
    f.write(f'{rows[7]}{i + 1},end_{barcode},{barcode}{end}\n')