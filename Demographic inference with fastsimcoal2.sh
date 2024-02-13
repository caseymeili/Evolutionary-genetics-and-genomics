#!/bin/bash

# Evolutionary Genetics and Genomics (HGEN 6092)
# R Session 4 - Population demography and structure (demographic inference with fastsimcoal2)
# Casey Meili

# -m 	using a MAF rather than DAF
# -n	number of coalescent simulations. Usually recommended to be 200,000-1,000,000
# -L	number of optimization (ECM) cycles to run to estimate the parameters Usually recommended to be between 50-100
# -M	indicates that parameter estimation will be performed
# -q	minimizes messages to the console
# -x	does not generate Arlequin output

./fsc28-mac -t BEI_expansion.tpl -e BEI_expansion.est -m -n 10000 -L 40 -M -q -x

2*3-2*(-201769.095/log10(exp(1)))
# substitute k with the # parameters (1 for constant, 3 for expansion, 4 for bottleneck)
# expansion model -> k-2 = 3-2
# -201769.095 = MaxEstLhood value from previous output 
