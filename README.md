## What is Botanist?

A quick analyzer for plant compounds, either the compounds themselves entered by name or id, or with a spectrometric input of a plant to break down (major) phytochemicals present with estimated percentages. It is 'quick' because it uses Woodward-Fieser rules in place of Density Functional Theory or other heavy computational methods, albeit compromising some accuracy. 

This particular tool computes an absorption or reflection spectrum for 120-800 nm, a bit broader than standard UV-Vis but just narrow enough to maintain fair accuracy for the methods used. 
Crucially, it relies on this dataset and is thus compound search is limited to what it contains: 
https://www.kaggle.com/datasets/yashasvigoswami/phytochemicals

Credits: Yashasvi Goswami (dataset owner), https://pubchem.ncbi.nlm.nih.gov/ (attributed data source)

# Capabilities 

analyze_compound: Takes Compound ID or compound name and plots absorption/reflection spectrum of a phytocompound. Predict plausible concentrations of phytocompounds from observed color of a plant part. Can download data. Using the respective flags, one may analyze a single or multiple compounds by id or by name. If multiple compounds are analyzed, they will 

analyze_plant: Opens a file selection window to choose the .csv and .pkl plant spectrometry files. Estimated percentage concentrations of phytocompounds in the plant part will be computed. Can download data.

# Unfinished

- cannot analyze multiple compounds yet
- failure to recognize spectral features for some compounds such as Acetic Acid and trace metals such as Gold

## Code Examples

[WIP]