# CCInclusive
CCInclusive studies anatree

## How to produce plots for the CC inclusive public note for Neutrino2016

On a Femilab gpvm machine (root6 is needed):

```
source /nusoft/app/alt/setup
setup root v6_04_06 -q e9:nu:prof
root -l HistoProducerModels.C
```

A file called histograms_TEM_MEC_trkrange_costheta_phi.root will be generated. This file contains the histograms to make the plots.
