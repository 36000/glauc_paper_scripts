import numpy as np
import altair as alt
import pandas as pd
from altair_transform import transform_chart
from altair import datum

pheno = pd.read_csv("output/pheno.csv", low_memory=False)

dis_subs = pheno[np.logical_not(
    np.isnan(pheno["6119-0.0"]) &
    np.isnan(pheno["6119-1.0"]) &
    np.isnan(pheno["6119-2.0"]) &
    np.isnan(pheno["6119-3.0"]))].eid.unique()

for acuity_f in ["5201-0.0", "5201-1.0", "5208-0.0", "5208-1.0"]:
    pheno = pheno[(pheno[acuity_f] <= 0.3) | np.isnan(pheno[acuity_f])]


pheno = pheno[
    (pheno["6148-2.0"] == -7)
    | (pheno["6148-3.0"] == -7)
    | pheno.eid.isin(dis_subs)]

profile = pd.read_csv("output/tract_profiles_wide.csv", low_memory=False)

spec_pheno = pd.DataFrame()
ethnic_map = {
    1: "White",
    1001: "British",
    2001: "White & Black Caribbean",
    3001: "Indian",
    4001: "Caribbean",
    2: "Mixed",
    1002: "Irish",
    2002: "White & Black African",
    3002: "Pakistani",
    4002: "African",
    3: "Asian or Asian British",
    1003: "Other white",
    2003: "White and Asian",
    3003: "Bangladeshi",
    4003: "Other Black",
    4: "Black or Black British",
    2004: "Other mixed",
    3004: "Other Asian",
    5: "Chinese",
    6: "Other ethnic group",
    -1: "Do not know",
    -3: "Prefer not to answer"}
sex_map = {
    0: "Female",
    1: "Male"}

pheno["21000-0.0"] = pheno["21000-0.0"].replace(np.nan, -1)
spec_pheno["Ethnicity"] = pheno["21000-0.0"].replace(ethnic_map)
spec_pheno["Age"] = pheno["21003-2.0"].astype(int)
spec_pheno["Sex"] = pheno["31-0.0"].replace(sex_map)
spec_pheno["TDI"] = pheno["189-0.0"].replace(np.nan, -2.13554).astype(int)
spec_pheno["subjectID"] = pheno["eid"]
spec_pheno = spec_pheno[spec_pheno["subjectID"].isin(pheno.eid)]
spec_pheno["is_dis"] = spec_pheno["subjectID"].isin(dis_subs)
spec_pheno["is_dis_disp"] = spec_pheno["is_dis"].map({True: 'Glaucoma', False: 'Control'})

charts = []
for ctrl in ["sec", "none"]:
    charts = []
    for ph_name, legend_loc in {"Age": [45, 81], "Sex": None, "Ethnicity": None, "TDI": None}.items():
        ctrl_subs = np.loadtxt(f"output/pheno_match/glauc_{ctrl}/ctrl_sub.txt")
        this_spec_pheno = spec_pheno[spec_pheno.is_dis | spec_pheno["subjectID"].isin(ctrl_subs)]

        total_count = len(this_spec_pheno.subjectID)
        print(f"Total # of subjects {ctrl}: {total_count}")

        for_violin_spec_pheno = spec_pheno.copy()
        if legend_loc is not None:
            for i in range(100):
                to_append = {
                    ph_name: 82,
                    "is_dis_disp": "100 subjects"
                }
                if i == 50:
                    to_append["texthere"] = True
                for_violin_spec_pheno = for_violin_spec_pheno.append(to_append, ignore_index=True)
            for i in range(200):
                for_violin_spec_pheno = for_violin_spec_pheno.append({
                    ph_name: 82,
                    "is_dis_disp": "White Space"
                }, ignore_index = True)

            y_spec = alt.Y(
                ph_name,
                scale=alt.Scale(domain=legend_loc))
        else:
            y_spec = alt.Y(
                ph_name)
        color_spec = alt.Color(
            'is_dis_disp:N', title="",
            scale=alt.Scale(
                domain=['Glaucoma', 'Control', "50 subjects", "100 subjects", "White Space"],
                range=['blue', 'orange', "black", "black", "white"]),
            legend=alt.Legend(values=["Glaucoma", "Control"]))
        x_spec = alt.X(
            f'count({ph_name}):Q',
            stack='center',
            impute=None,
            title=None,
            axis=alt.Axis(labels=False, values=[0], grid=False, ticks=True)
        )

        if legend_loc is None:
            title = ""
            fontsize = 40
        else:
            title = "A. Unmatched"
            fontsize = 20
        chart2 = alt.Chart(for_violin_spec_pheno, title=title).mark_bar(orient='horizontal').encode(
            y=y_spec,
            color=color_spec,
            x=x_spec
        )

        chart_text = alt.Chart(for_violin_spec_pheno).mark_text(
            align='left',
            baseline='top',
            fontSize = fontsize,
            dx = 18,
            dy = 5
        ).encode(
            x=x_spec,
            y=y_spec,
            text='is_dis_disp',
            color=color_spec
        ).transform_filter(datum.texthere)

        chart2 = chart2 + chart_text

        for_violin_spec_pheno = this_spec_pheno.copy()
        if legend_loc is not None:
            for i in range(50):
                to_append = {
                    ph_name: 82,
                    "is_dis_disp": "50 subjects"
                }
                if i == 25:
                    to_append["texthere"] = True
                for_violin_spec_pheno = for_violin_spec_pheno.append(to_append, ignore_index=True)
            for i in range(100):
                for_violin_spec_pheno = for_violin_spec_pheno.append({
                    ph_name: 82,
                    "is_dis_disp": "White Space"
                }, ignore_index = True)

        if legend_loc is None:
            title = ""
        else:
            title = "B. Matched"
        chart1 = alt.Chart(for_violin_spec_pheno, title=title).mark_bar(orient='horizontal').encode(
            y=y_spec,
            color=color_spec,
            x=x_spec
        )

        chart_text = alt.Chart(for_violin_spec_pheno).mark_text(
            align='left',
            baseline='top',
            fontSize = fontsize,
            dx = 18,
            dy = 5
        ).encode(
            x=x_spec,
            y=y_spec,
            text='is_dis_disp',
            color=color_spec
        ).transform_filter(datum.texthere)

        chart1 = chart1 + chart_text

        if not ctrl == "none":
            chart1 = chart2 | chart1
        if legend_loc is None:
            charts.append(chart1)
        else:
            chart1.save(f'output/second_paper_plots/pheno_hist_{ctrl}_{ph_name}.html')
    alt.VConcatChart(vconcat=charts).save(f'output/second_paper_plots/pheno_hist_{ctrl}.html')

