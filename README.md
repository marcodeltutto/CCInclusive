# CCInclusive
CCInclusive studies anatree

## How to produce plots for the CC inclusive public note for Neutrino2016

### Model variations

On a Femilab gpvm machine (root6 is needed):

```
source /nusoft/app/alt/setup
setup root v6_04_06 -q e9:nu:prof
root -l HistoProducerModels.C        # for MEC && ESF+TEM
root -l HistoProducerMA.C            # for Ma CCQE
```

A file called `histograms_TEM_MEC_trkrange_costheta_phi.root` will be generated for MEC && ESF+TEM and one called `histograms_MA_trkrange_costheta_phi.root` for Ma CCQE. This files contain the histograms to make the plots.

### How to calculate the simulated POT

Open `calculatePOT.py` and change the name of the file, then run:
```
python calculatePOT.py
```
