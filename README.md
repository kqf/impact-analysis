# Impact analysis ![status](https://travis-ci.com/kqf/impact-analysis.svg?token=7bkqqhrPB19pD1YKrAZM&branch=master)


This repo contains the code that estimates the profile of the impact–parameter dependent elastic scattering amplitude. In original analysis the input dataset included the experimental data on elastic hardron-hadron scattering starting from 28 GeV (ISR) up to the LHC energies (7 TeV). Impact-parameter dependent inelastic overlap function helps to reveal the asymptotics of hadron interactions. In particular LHC data should improve our understanding of the mechanisms at high energies.

## Usage
For some details/examples check the [tests](test) directory.
```bash

# First step:
# 			Make sure that you have the valid dataset in the input/*.data

# Second step:
#			Select the energies of your interest

vim config/input.json

# Third step: 
#			Run the analysis

make

```
## NB
Original version of the code is lost (laptop damage + didn't push the code). This version of the draft works well and produces the same results as in the paper.
