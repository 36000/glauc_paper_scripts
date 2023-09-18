
import pandas as pd
import numpy as np
import altair as alt
import scipy.stats as st

datasets = ["glauc_sec", "age", "amd_rob", "amd_rob_age", "glauc_rob", "tenano"]
bundles = ["UNC_L", "UNC_R", "CST_L", "CST_R", "L_OR", "R_OR"] #"fov_L", "fov_R",  "mac_L", "mac_R", "periph_L", "periph_R"]
formal_b_names = {
    "UNC_L": "UNC L",
    "UNC_R": "UNC R",
    "CST_L": "CST L",
    "CST_R": "CST R",
    "L_OR" : "OR L",
    "R_OR" : "OR R",
    "fov_L": "fOR L",
    "fov_R": "fOR R",
    "mac_L": "mOR L",
    "mac_R": "mOR R",
    "periph_L": "pOR L",
    "periph_R": "pOR R",
}

prof_csvs = pd.DataFrame()
for ds in datasets:
    this_df = pd.read_csv(f"output/pheno_match/{ds}/profiles.csv")
    this_df = this_df[this_df["nodeID"] == 49]
    this_df["dataset_name"] = ds
    prof_csvs = pd.concat([prof_csvs, this_df], ignore_index=True)

prof_csvs["dataset_name"] = prof_csvs["dataset_name"].replace(
    {
        "glauc_sec": "A",
        "age": "B",
        "amd_rob": "A.1",
        "amd_rob_age": "B.2",
        "glauc_rob": "B.1",
        "tenano": "A.2"
    }
)
reduced_prof_csvs = {}
reduced_prof_csvs["dis_status"] = [0]*len(prof_csvs) + [1]*len(prof_csvs)
reduced_prof_csvs["eid"] = list(prof_csvs["subjectID_ctrl"]) + list(prof_csvs["subjectID_dis"])
reduced_prof_csvs["dataset_name"] = list(prof_csvs["dataset_name"]) + list(prof_csvs["dataset_name"])
for bundle in bundles:
    reduced_prof_csvs["ep"+bundle] = list(prof_csvs[f"dki_fa_{bundle}_ctrl"]) + list(prof_csvs[f"dki_fa_{bundle}_dis"])
reduced_prof_csvs = pd.DataFrame(reduced_prof_csvs)

for bundle in bundles:
    reduced_prof_csvs["ep"+bundle] = ~reduced_prof_csvs["ep"+bundle].isna()
    print(f"{bundle}: {np.sum(reduced_prof_csvs['ep'+bundle])/len(reduced_prof_csvs['ep'+bundle])}")
reduced_prof_csvs = pd.wide_to_long(
    reduced_prof_csvs, "ep", ["dis_status", "dataset_name", "eid"], "bname", suffix="\D+").reset_index()
reduced_prof_csvs["dis_status"] = reduced_prof_csvs["dis_status"].replace({0: ", Control", 1: ", Test"})
reduced_prof_csvs["bname"] = reduced_prof_csvs["bname"].replace(formal_b_names)
reduced_prof_csvs["Bundle, label"] = reduced_prof_csvs["bname"] + reduced_prof_csvs["dis_status"]

reduced_prof_csvs["Percent Found"] = reduced_prof_csvs["ep"]
reduced_prof_csvs["Dataset Name"] = reduced_prof_csvs["dataset_name"]

chart = alt.Chart().mark_bar().encode(
   color=alt.Color("Bundle, label", scale=alt.Scale(scheme="category20")),
   y=alt.Y("mean(Percent Found):Q"))

error_bars = alt.Chart().mark_errorbar(extent='ci').encode(
    y=alt.Y('Percent Found:Q')
)

font_size = 40
alt.LayerChart(
    layer=[chart, error_bars],
    data=reduced_prof_csvs).facet(
        column=alt.Column("Bundle, label", header=alt.Header(labelAngle=90)),
        row=alt.Row("Dataset Name", header=alt.Header(labelFontSize=font_size))).configure_title(
        fontSize=font_size).configure_axis(
            labelFontSize=20,
            titleFontSize=font_size,
            labelLimit=0).configure_legend(
                labelFontSize=font_size,
                titleFontSize=font_size,
                labelLimit=0,
                titleLimit=0,
                columns=2,
                orient='right'
            ).resolve_scale(y='shared').save('output/second_paper_plots/missing_bias.html')
