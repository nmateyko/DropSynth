'''
Calculates number of oligos per bead given bead diameter in μm and oligo concentration in μM.
'''
import sys
from decimal import Decimal

diameter = float(sys.argv[1])
concentration = float(sys.argv[2])

volume_ul = 1e9 * ((4 / 3) * ((((diameter / 2) / 1e6) ** 3) * 3.14))
oligo_per_bead = 6.02e23 * concentration * volume_ul * 1e-12

print('%.2E' % Decimal(oligo_per_bead))