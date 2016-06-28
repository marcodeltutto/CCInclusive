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
Or you can avoid re-generating them and use the ones in this repository.

### How to calculate the simulated POT

Open `calculatePOT.py` and change the name of the file, then run:
```
python calculatePOT.py
```

### How to produce the plots

Run:
```
root -l "draw_MA.C(selection,variable)"
root -l "draw_TEM_MEC.C(selection,variable)"
```
where selection = 1 is Christoph selection and selection = 2 is Xiao selection. Variable is the quantity you want on the x axis: 0 = track length, 1 = cos(theta), 2 = phi.
It may require the rootlogon.C file in `/nashome/m/mdeltutt/rootlogon.C`.

