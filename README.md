# A method to detect a sudden increase in prevelance after introduction of outpatient registry

### Objective

To detect, for each disease endpoint, whether there is a significant increase in prevalence in 1998 likely due to introduction of outpatient register

![alt text](https://github.com/dsgelab/jump1998/image/objective.png)

- Because many disease outcomes are missed before 1998:
 - Analyses that consider age at onset are not trustable
 - Bias is introduced due to swift change in prevalence

- Dataset: 
 - Nationwide Finnish health registry data (~ 5.7 Million individuals)
 - Construct the same definitions of the FinnGen endpoints

### We use Poisson Regression to model changes in  disease prevalence

For each disease - 
- Covariates:
 - event year
 - age
 - “discontinuation indicator”:
    1 if the event year falls into 1999-2001
    0 if the event year falls into 1995-1997
- Output: prevalence of the disease
- Time frame: 1986 - 2011
