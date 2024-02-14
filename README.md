First generate the tract profiles by following the steps outlined here:
    https://github.com/36000/OR_aging_ukbb

To process individual subjects, follow all of the instructions in the
afq_processing folder of that repository.

To generate the tract profiles, follow the steps in the 
analysis_and_figures folder up to running gen_pheno.py .

You can also use the python_analysis_docker from that repository to
run these scripts.

For some of these plots, you can also skip processing of all of the subjects, and instead use the aggregrate data provided in `aggregate_data`. For others, the models that generate the results are provided in `ML_models`.

For figure 1, run: first_paper_plots/general_profiles_plot.py glauc_match_profs_left

**_NOTE:_** The figure 1 script can also be run on `profiles_w_aci_glauc_sec_match.csv` in `aggregate_data`.

For figures 2 and 3, run primary_match.py with `pheno_of_interest` set to `"glauc_sec"`

**_NOTE:_** The models used to generate figures 2 and 3 are provided in `ML_models` under `glaucoma`.

For figure 4, run these series of commands:
  * primary_match.py with `pheno_of_interest` set to `"age"`
  * secondary_match.py amd_rob
  * secondary_match.py tenano
  * secondary_match.py glauc_rob
  * secondary_match.py amd_rob_age

**_NOTE:_** The models used to generate figure 4 are provided in `ML_models` under `age`.

For supplementary figures 1 and 2, run pheno_histograms.py

**_NOTE:_** Supplementary figures 1 and 2 require access to UK Biobank phenotypic data which cannot be made publicly available.

For supplementary figure 3, run first_paper_plots/general_profiles_plot.py glauc_match_profs_right

**_NOTE:_** This script can also be run on `profiles_w_aci_glauc_sec_match.csv` in aggregate data.

For supplementary figure 4, run missing_bias.py

**_NOTE:_** Supplementary figure 4 requires access to individual tract profiles, which cannot be made publicly available.

For supplementary figure 5, run full_bundle_test.py

**_NOTE:_** Supplementary figure 5 uses models provided in `ML_models` under `glaucoma`.

Notice that there are many options for `pheno_of_interest` in primary_match.py.
This includes the options `"glauc_sec_no_cut"` and `"age_nocut"`, which test
whether the main results of the paper hold when not cutting subjects with
high logmar values. It also includes options for doing statistical
matching with other control variables, other cutoff values, and other
techniques. In all cases, the main results of the paper hold.
