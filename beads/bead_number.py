'''
Calculates number of beads per mL given bead diameter in μm and packing efficiency in percent.
Default packing efficiency is 63.5% (https://en.wikipedia.org/wiki/Sphere_packing)
'''
import sys
from decimal import Decimal

diameter = int(sys.argv[1])
if len(sys.argv) > 2:
  efficiency = int(sys.argv[2])
else:
  efficiency = 63.5

bead_per_mL = (1e12 * (efficiency / 100)) / ((4 / 3) * ((diameter / 2) ** 3) * 3.14)

print('%.2E' % Decimal(bead_per_mL))