#!/usr/bin/python
import sys
import simulate_prewetting

c_bound = float(sys.argv[1])
Js = float(sys.argv[2])
chem_potent1 = float(sys.argv[3])
chem_potent2 = float(sys.argv[4])
fidx = int(sys.argv[5])
speedup=float(sys.argv[6])
simulate_prewetting.main(c_bound, Js, chem_potent1, chem_potent2, fidx,speedup)
