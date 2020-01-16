# LPM-Geomod_Suite
This is a geologic modeling suite plugins for AR2GEMS using signed distances functions and indicators.

## Installation
* Run *install_dependencies.py* script on python terminal in AR2GEMS folder;  
* Paste the ui folder content to *plugins\Geostat* AR2GEMS folder; 
* Paste py folder content to *plugins\Geostat\python* AR2GEMS folder; 
* Paste *helpers.py* in Lib AR2GEMS folder.

If everything goes well you will see LPM-Geomod_Suite in the algorithms pannel.

![Algo pannel](images/algo_pannel.PNG)

## Transformations 
Tranformations are straightforward. User needs to select the categorical rock type property.

### Indicators
Creates one indicator property for each category.

### Signed distances
Creates one signed distances property for each property.

## Deterministic
This plugin will create a deterministic geological model from indicators or signed distances properties. The interpolation method is dual kriging which is global it can take too long on big datasets. iterations options will refine the grid defining a refinement zone from a marching cube.

### Transformations
There are transformations for the gridded properties too.

#### Softmax transfomration

Softmax transforms interpolated distances into probabilities or correct order relations issues with interpolated indicators properties. The gamma parameter control the magnitude of the uncertainty.

#### Entropy

Entropy is calculated from probabilities properties. It returns a single property that measures multicategorical uncertainty. It basically delineates the contacts.

## Stochastic

## Validation

## Optimizationn

## VTK export

This plugin exports point sets and/or grids in VTK format. Esported files can be opened with [Paraview](http://www.paraview.org).

## Auto grid

Auto grid creates automatically a grid which covers the point set given the block dimensions and a buffer size.

## What can be done?
