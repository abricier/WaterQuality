# WaterQuality

[![DOI](https://zenodo.org/badge/224832878.svg)](https://zenodo.org/badge/latestdoi/224832878)

## Synopsis

The WaterQuality package extends the functionalities of the WaterDetect package to calculate continental water quality parameters from satellite reflectances. It has been specially conceived for L2A Sentinel 2 imagery from [MAJA](https://logiciels.cnes.fr/en/content/maja)<sup>1</sup>  processor, and the parameter is calculated just where there exist water, according to the waterdetect mask. Several inversion algorithms from the literature have been implemented:<br>
* Chlorophyll - Lins
* Chlorophyll - Gitelson
* CDOM absorption - Brezonik
* Turbidity - Dogliotti
* SPM - adapted from Nechad at the GET laboratory

![Screenshot](GraphicalAbstract.JPG)

All the details and tests has been described in the article <b>Automatic Water Detection from Multidimensional Hierarchical Clustering for Sentinel-2 Images and a Comparison with Level 2A Processors</b>, under revision by the journal Remote Sensing of Environment.

<b>How to cite ("accepted by Remote Sensing of Environment, pending publication"):</b><br>
Cordeiro, M.C.R, Martinez, J.-M., Pena Luque, S., 2020. Automatic Water Detection from Multidimensional Hierarchical Clustering for Sentinel-2 Images and a Comparison with Level 2A Processors. Remote Sensing of Environment XX, XX. 


## Dependencies
The required libraries are:
```
GDAL>=3.0.2
matplotlib>=3.1.2
PyPDF2>=1.26.0
scipy>=1.3.2
scikit-learn>=0.22
skimage>=0.16.2
numpy>=1.17
waterdetect>=1.5
```
The test_dependencies.py can be used to check if all libraries are loading correctly. Simply run:

```python test_dependencies.py```

## Usage
To use it, you should clone the project to your repository and run "python runWaterColor.py --help"
```
usage: runWaterColor.py [-h] -i INPUT -o OUT [-s SHP] [-p PRODUCT] [-g]
                        [-c CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The products input folder. Required.
  -o OUT, --out OUT     Output directory. Required.
  -s SHP, --shp SHP     SHP file. Optional.
  -p PRODUCT, --product PRODUCT
                        The product to be processed (S2_Theia, Landsat,
                        S2_L1C)
  -g, --off_graphs      Turns off the scatter plot graphs
  -c CONFIG, --config CONFIG
                        Configuration .ini file. If not specified
                        WaterDetect.ini is used as default
```

The input directory should contain the uncompressed folders for the images. The script will loop through all folders in the input directory and save the water masks, graphs and reports to the output folder. The output folder must be created beforehand.

If the config file is not specified, the script will search for WaterDetect.ini in the current folder.

## Config File

## Contributors
> Author: Maurício Cordeiro (ANA/GET)<br>
> Supervisor: Jean-Michel Martinez (IRD/GET)<br>
> Validation dataset: Santiago Pena Luque (CNES) 

### Institutions
* ANA - Agência Nacional de Águas (https://www.gov.br/ana/en/)
* GET - Géosciences Environnement Toulouse (https://www.get.omp.eu/)
* IRD - Institut de Recherche pour le Développement (https://en.ird.fr/)
* CNES - Centre National d'Études Spatiales (https://cnes.fr/fr)

## License
This code is licensed under the [GNU General Public License v3.0](https://github.com/cordmaur/WaterDetect/blob/master/LICENSE) license. Please, refer to GNU's webpage  (https://www.gnu.org/licenses/gpl-3.0.en.html) for details.

## Reference
(1) Hagolle, O.; Huc, M.; Pascual, D. V.; Dedieu, G. A Multi-Temporal Method for Cloud Detection, Applied to FORMOSAT-2, VENµS, LANDSAT and SENTINEL-2 Images. Remote Sensing of Environment 2010, 114 (8), 1747–1755. https://doi.org/10.1016/j.rse.2010.03.002.
