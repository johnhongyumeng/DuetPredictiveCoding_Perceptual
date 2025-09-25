# DuetPredictiveCoding_Perceptual

Code for: Duet model unifies diverse neuroscience experimental findings on predictive coding
John H. Meng, Jordan M. Ross, Jordan P. Hamm, Xiao-Jing Wang
BioRxiv link: https://doi.org/10.1101/2025.07.12.664417

The code is published under MIT license.
To run, clone or download the whole repository. The code is tested under Python 3.9 and Matlab R2024b.

For a quick demo, just run
./scripts/class_main_oddball.py   running time： ~1 min


## Figure 2: Omission
./scripts/Shell_main_interval_omission.py 
Requires: class_main_oddball

Please note that due to size limitations, to run the code with the realistic version (i.e., without homeostasis assumption, it is needed to run the code from https://github.com/johnhongyumeng/DuetPredictiveCoding_Plasticity 
to initiate the local connectivity. 

## Figure 3:  Matching to MMN
./scripts/class_main_oddball_prob.py.    
Running time: 7  min
## Figure 4: Number of repetitions.
./scripts/Shell_main_nRep.py  
Which requires
class_main_oddball_prob_nRep_dev.py 10 mins

## Figures 5: Fixed oddball 
Run: ./scripts/class_main_oddball.py
## Figure 6: Randomized oddball
Run: ./scripts/class_main_oddball_prob.py
## Figure 7: Analyzing code
Run: ./analyze/



Data analysis:
To run the code, the original data is needed from Bastos et al., (2023), which can be accessed through

https://zenodo.org/record/8335415

and

https://doi.org/10.5281/zenodo.8335415/ 

After accessing the data, run

./analyze/BATCH_analyse_modified.m 

./analyze/BATCH_analyse_plot.m 

for normalized traces.

and run

./analyze/Analyse_meng_tuningcurve_combine_BH_tuning.m 

for tuning curves.

Please remember to check the commented lines for proper modifications. 
