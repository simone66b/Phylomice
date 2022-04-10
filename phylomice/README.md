# phylomice
Extension of R [mice](https://github.com/stefvanbuuren/mice) package adding new imputation methods that incorporate phylogenetic information.

You should be familiar with the mice package to use this one.

## Installation
You can install the package using [devtools](https://github.com/hadley/devtools):
```
devtools::install_github("pdrhlik/phylomice")
```
### phnorm
Imputes univariate continuous missing data using the generalized least square approach. It is based on the **mice.impute.norm** method in **mice** package.
```
library(phylomice)
imp <- mice(data, method = 'phnorm', tree=tree)
```
### phpmm
Uses predictive mean matching from the pGLS regression line, and biases selection towards more closely related species.

## Authors

* **Patrik Drhlik** - *Initial work* - [pdrhlik](https://github.com/pdrhlik)
* **Simone P. Blomberg** - *Theoretical background*
 
If you are interested in improving current or creating new imputation methods that use phylogenies, don't hesitate to contribute.

## License
This project is licensed under the GPL-3 License.

## Acknowledgments
This project wouldn't be possible without the following parties:
* The University of Queensland, Australia
* Technical University of Liberec, Czech Republic
* NESSIE Erasmus Mundus project
