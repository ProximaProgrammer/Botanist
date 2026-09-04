# What is Botanist?

A quick analyzer for plant compounds, either the compounds themselves entered by name or id, or with a spectrometric input of a plant to break down (major) phytochemicals present with estimated percentages. It is 'quick' because it uses Woodward-Fieser rules in place of Density Functional Theory or other heavy computational methods, albeit compromising some accuracy. Component identification algorithms were inspired by the dataset used in https://zenodo.org/records/17520032

This particular tool computes an absorption or reflection spectrum for 120-800 nm, a bit broader than standard UV-Vis but just narrow enough to maintain fair accuracy for the methods used. 
Crucially, it relies on the following dataset and is thus compound search is limited to what it contains: \
\
https://www.kaggle.com/datasets/yashasvigoswami/phytochemicals \
(dataset owner: Yashasvi Goswami, data source: https://pubchem.ncbi.nlm.nih.gov/)

## Capabilities 

***analyze_compound***: Takes Compound ID or compound name and plots absorption/reflection spectrum of a phytocompound. Predict plausible concentrations of phytocompounds from observed color of a plant part. Can download data. Using the respective flags, one may analyze a single or multiple compounds by id or by name. If multiple compounds are analyzed, they will 

***analyze_plant***: Opens a file selection window to choose the .csv and .pkl plant spectrometry files. Estimated percentage concentrations of phytocompounds in the plant part will be computed. Can download data.

## Unfinished

- cannot analyze multiple compounds yet
- only absorption spectrum, no reflection spectrum
- limited 

- failure to recognize spectral features for some compounds such as Acetic Acid and trace metals such as Gold (though this is fundamentally a cost using the 'shortcut' methods this project is founded on)

# Setup

```
$ git clone https://github.com/ProximaProgrammer/Botanist
$ cd Botanist

$ mkdir build && cd build
$ cmake ..
$ cmake --build . --config Release
$ cd ..

$ pip install .
```

# Code Examples

`$ botanist --help`          
```
usage: botanist [-h] [-analyze_compound] [-hide] [-d {csv,txt,json}] [-n NAME]
                [--names NAMES [NAMES ...]] [-id CID] [-ids CIDs [CIDs ...]]
                [-spectrum_type {absorption,reflection}] [--analyze_plant]
                ...
```
\
`$ botanist -analyze_compound -n Caffeine -spectrum_type absorption`
```
Electrons: total  102  | valence  74  | core  28
=== ESTIMATES ===
ketone_enone | lambda=272 nm | double_bonds=4 | atoms=11
  -> spectral points: 242
--displaying--
[this should display a graph for the single compound specified]
```
\
`$ botanist -analyze_compound -n "(-)-10-EPI-ALPHA-CYPERONE" -spectrum_type absorption`
```
Electrons: total  120  | valence  88  | core  32
=== ESTIMATES ===
ketone_enone | lambda=242 nm | double_bonds=2 | atoms=4
  -> spectral points: 127
--displaying--
[this should display a graph for the single compound specified]
```
\
`$ botanist -analyze_plant`
```
--displaying--
{'Total_Chlorophyll': 0.159, 'Total_Carotenoids': 0.012, 'Anthocyanins': 0.0, 'Total_Phenolics': 0.0}
```


![WIP](sussy.gif)

> Built by Prox (Discord: wprox)