Code and data for the paper [The functional role of sequentially neuromodulated synaptic plasticity in behavioural learning](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009017)

## [Experimental data](data/experiments)

Variables in `trial_data.csv` : 
* `sub`: subject ID 
* `phase`: task stage. 0 for initial learning, 1 for initial learning 
* `condition`: experimental group. 'OFF' - light-off control, 'ON' - light-on ACh-suppressed group, 'GFP' - light-on control 
* `trial`: recoded trial variable. 1.1 for trial 1 on day 1, 1.2 for trial 2 on day 1 ...
* `rew.found`: 1 if found, 0 if not 
* `day`: 1 - 8 for initial learning, 1 - 12 for reversal learning 

Processed and summarized data can be accessed in `processed_trial_data.rds` using R

```
> datalist <- readRDS("../../data/experiments/processed_trial_data.rds")
> names(datalist)
[1] "combined.dat"    "mouse.success.p" "group.success.p" "ids"            
[5] "hit.80.ids"      "hit.80"
```
* `combined.dat`: raw data grouped by subject ID in a tibble
* `mouse.success.p`: successful trials per day of each mouse
* `group.success.p`: group average performance 
* `hit.80.ids`: number of days taken to reach 80% performance with subject IDs

## Code 

### [sn-Plast model](code/snPlast_model) 

Under `/snPlast_model`:
* Python implementation of the [sn-Plast model in Brzosko et al., 2017](https://elifesciences.org/articles/27756) adapted to the experimental open-field task found in the `snPlast_model` folder. 
* `param_search.sh`: example bash script to run the model across different parameter settings (e.g. <img src="https://latex.codecogs.com/gif.latex?\eta_{ACh}\text{,}\eta_{DA}" /> )
