First generate the tract profiles by following the steps outlined here:
    https://github.com/36000/OR_aging_ukbb

To process individual subjects, follow all of the instructions in the
afq_processing folder of that repository.

To generate the tract profiles, follow the steps in the 
analysis_and_figures folder up to running gen_pheno.py .

You can also use the python_analysis_docker from that repository to
run these scripts.

For figure 1, run: first_paper_plots/general_profiles_plot.py glauc_match_profs_left

For figure 2, run primary_match.py with `pheno_of_interest` set to `"glauc_sec"`

For figure 3, run primary_match.py with `pheno_of_interest` set to `"age"`

For figure 4, run these series of commands:
  * secondary_match.py amd_rob
  * secondary_match.py tenano
  * secondary_match.py glauc_rob
  *  secondary_match.py amd_rob_age

For figure 5, run shap_comp.py

For supplementary figures 1 and 2, run pheno_histograms.py

For supplementary figure 3, run missing_bias.py

For supplementary figure 4, run full_bundle_test.py

Notice that there are many options for `pheno_of_interest` in primary_match.py.
This includes the options `"glauc_sec_no_cut"` and `"age_nocut"`, which test
whether the main results of the paper hold when not cutting subjects with
high logmar values. It also includes options for doing statistical
matching with other control variables, other cutoff values, and other
techniques. In all cases, the main results of the paper hold.
