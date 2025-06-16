### Modeling codon usage and tRNA dynamics

This repo contains code and data associated with the following paper: A. M. Hill, K. To, C. O. Wilke (preprint). Availability of charged tRNAs drives maximal protein synthesis at intermediate levels of codon usage bias.

The first part of this paper considers a mathematical model for tRNA charging dynamics than can be solved numerically. The analysis can be reproduced by running each Python notebook in `mathematical_analysis` (requires pandas, seaborn, and sympy)

The second part is a simulation study which was done using the [Pinetree](https://github.com/clauswilke/pinetree) stochastic simulation software (version 0.4.1). More info about pinetree can be found here: https://pinetree.readthedocs.io/en/latest/

To install pinetree, download the source code from https://github.com/clauswilke/pinetree. The easiest way to install is using pip:

```
pip3 install cmake   
git clone https://github.com/clauswilke/pinetree.git
cd pinetree
pip3 install .
```

Python scripts for the simulations are located in `src/python/models`. There are two models---`revised_weighted_phage_model.py` sets up a simulation of bacteriophage T7 infection with static tRNAs/codon speeds, and `trna_phage_model.py` sets up an identical simulation, but with dynamic tRNAs. These can both be run from the command line (with pinetree installed):

```
python3 trna_phage_model.py <fop> <charge_rate> <pref_proportion> <seed_val> <ribo_speed> <trna_count>
```

- `fop` is the fraction of optimal codons in T7 gene 10A (should be a number from 0 to 1)
- `charge_rate` specifies the reaction rate constant for tRNA recharging
- `pref_proportion` is the proportion of tRNAs that are preferred (which correspond to optimal codons). This is also a number from 0 to 1.
- `seed_val` seed to reproduce simulations
- `ribo_speed` the baseline elongation rate for ribosomes in the simulation
- `trna_count` the total number of tRNAs in the simulation

```
python3 revised_weighted_phage_model.py <fop> <speed_op> <speed_non_op> <seed_val> <ribo_speed>
```

- `fop` is the fraction of optimal codons in T7 gene 10A (should be a number from 0 to 1)
- `speed_op` static codon speed associated with optimal codons
- `speed_non_op` static codon speed associated with non-optimal codons
- `seed_val` seed to reproduce simulations
- `ribo_speed` the baseline elongation rate for ribosomes in the simulation

Finally, code to produce all figures in the manuscript is in `R`.
