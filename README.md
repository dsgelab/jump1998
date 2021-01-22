## A method to detect a sudden increase in prevelance after introduction of outpatient registry

### Objective

To detect, for each disease endpoint, whether there is a significant increase in prevalence in 1998 likely due to introduction of outpatient register

![alt text](https://raw.githubusercontent.com/dsgelab/jump1998/main/images/objective.png)

- Because many disease outcomes are missed before 1998:
	* Analyses that consider age at onset are not trustable
	* Bias is introduced due to swift change in prevalence

- Dataset: 
	* Nationwide Finnish health registry data (~ 5.7 Million individuals)
	* Construct the same definitions of the FinnGen endpoints

### We use Poisson Regression to model changes in  disease prevalence

For each disease, 
- Covariates:
	* event year
	* age
	* “discontinuation indicator”: 1 if the event year falls into 1999-2001; 0 if the event year falls into 1995-1997
- Output: prevalence of the disease
- Time frame: 1986 - 2011

![alt text](https://raw.githubusercontent.com/dsgelab/jump1998/main/images/model.png)

### Experimental setting to test model performances

![alt text](https://raw.githubusercontent.com/dsgelab/jump1998/main/images/experiment1.png)

Our model accurately identify endpoints with significant prevalence increase in 1998

![alt text](https://raw.githubusercontent.com/dsgelab/jump1998/main/images/experiment2.png)

Percentage of endpoints having a significant prevalence increase in 1998 is 70.6%

![alt text](https://raw.githubusercontent.com/dsgelab/jump1998/main/images/experiment3.png)

Percentage of endpoints with a significant prevalence increase in Neoplasms - benign neoplasms vs malignant neoplasms: 60.1% vs 2.4%.
![alt text](https://raw.githubusercontent.com/dsgelab/jump1998/main/images/experiment4.png)
Reassuring results because Malignant neoplasm are well covered by the cancer register since the 70’s.


### A warning when a substantial prevalence change in 1998 is detected has been added to Risteys.
https://risteys.finngen.fi
An example: 
![alt text](https://raw.githubusercontent.com/dsgelab/jump1998/main/images/risteys.png)
