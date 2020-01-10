# LPM-Geomod_Suite
This is a geologic modeling suite plugins for AR2GEMS using signed distances functions and indicators.

## Installation
Run install_dependencies.py script from run script inside AR2GEMS.  
Paste the ui folder content to plugins\Geostat AR2GEMS folder and py folder content to plugins\Geostat\python AR2GEMS folder. helpers.py should go on Lib AR2GEMS folder.
If everything goes well you will see LPM-Geomod_Suite in the algorithms pannel.

![Algo pannel](images/algo_pannel.png)

## Transformations 
Tranformations are straightforward. User needs to select the categorical rock type property.
![Transformations](images/transform.PNG)

### Indicators
Creates one indicator property for each category.

### Signed distances
Creates one signed distances property for each property.

## Deterministic
This plugin will create a deterministic geological model from indicators or signed distances properties. The interpolation method is dual kriging which is global it can take too long on big datasets. iterations options will refine the grid defining a refinement zone from a marching cube.

User must choose the variables type (indicators or signed distances), add and select a variogram model and a code for each category.

### Transformations
There are transformations for the gridded properties too.

#### Softmax

#### Entropy

## Stochastic

## Validation

## Optimizationn

## VTK export

## What can be done?
