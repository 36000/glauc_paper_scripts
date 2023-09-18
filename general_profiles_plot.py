import numpy as np
import altair as alt
import pandas as pd
import os.path as op
import sys
from tqdm import tqdm

from global_configs import *

dataset = sys.argv[1] # can be hcp_roi, hcp_reco, , ukbb_r, ukbb_lvr, hcp_roi_*, hcp_reco_*

# dataset 1: ukbb_l_filt_oab

if "lr_only" in dataset:
    lr_only = True
else:
    lr_only = False

b_name = ""
f_pre = ""
use_iqr = False
if "sex" in dataset:
    binning = "sex"
elif "acuity_l" in dataset:
    binning = "acuity_l"
elif "site" in dataset:
    binning = "site"
else:
    binning = None

if "filt" in dataset:
    filt = True
else:
    filt = False

if "oab" in dataset:
    oab = True
else:
    oab = False

if "ukbb" in dataset:
    if "lvr" in dataset:
        sides = ["LvR"]
    elif "r" in dataset:
        sides = ["R"]
    elif "l" in dataset:
        sides = ["L"]
    else:
        raise ValueError("improper UKBB name")
    all_sides = False
    if sides[0] == "LvR":
        if "cst" in dataset:
            columns = ["CST_L", "CST_R", "LR_CST"]
            b_name = "_cst"
        elif "unc" in dataset:
            columns = ["UNC_L", "UNC_R", "LR_UNC"]
            b_name = "_unc"
        else:
            columns = ["LR_fov", "LR_mac", "LR_per"]
        no_LR = False
    else:
        columns = ["fov", "mac", "periph", "crossfp", "crossfm"]
        no_LR = True
    dataset = "ukbb"
elif "match" in dataset:
    no_LR = False

    if "glauc" in dataset:
        f_pre = "glauc_sec_"
    if "ambly" in dataset:
        f_pre = "ambly_"
    if "logmar" in dataset:
        f_pre = "logmar_"
    if "any" in dataset:
        f_pre = "any_"

    if "profs" in dataset:
        all_sides = False
        use_iqr = True
        binning = "is_dis"
        if "right" in dataset:
            columns = ["R_OR", "CST_R", "UNC_R"]
            sides = ["R"]
        else:
            columns = ["L_OR", "CST_L", "UNC_L"]
            sides = ["L"]
        dataset = "profs"
    else:
        all_sides = True
        columns = ["OR", "CST", "UNC"]
        sides = ["L", "R"]
        dataset = "match"
else:
    all_sides = False
    if dataset[-3:] == "lvr":
        sides = ["LvR"]
    elif dataset[-1] == "r":
        sides = ["R"]
    elif dataset[-1] == "l":
        sides = ["L"]
    else:
        all_sides = True
    if "hcp_reco" in dataset:
        dataset = "hcp_reco"
    elif "hcp_roi" in dataset:
        dataset = "hcp_roi"
    elif "hcp_bootroi" in dataset:
        dataset = "hcp_bootroi"
    else:
        dataset = "hcp_boot"
    if lr_only:
        columns = ["LR_only", "LR_cross"]
        sides = ["L", "R"]
        no_LR = False
    else:
        if all_sides:
            columns = ["fov", "periph", "LR", "cross"]
            sides = ["L", "R"]
            no_LR = False
        else:
            if sides[0] == "LvR":
                columns = ["LR_fov", "LR_per"]
                no_LR = False
            else:
                columns = ["fov", "periph", "cross"]
                no_LR = True

if binning is not None and binning != "is_dis":
    f_suf = f"_{binning}_binned"
else:
    f_suf = ""

if filt:
    f_suf = f_suf + "_filt"
if oab:
    f_suf = f_suf + "_oab"

profiles = pd.read_csv(f"output/aci_for_paper/profiles_w_aci_{f_pre}{dataset}{f_suf}.csv")
profiles.dropna(inplace=True)

if dataset != "match":
    for column_name in profiles:
        if "aci" not in column_name and "md" in column_name:
            profiles[column_name] = profiles[column_name] * 1000
profiles["0_LINE"] = 0
profiles["Left"] = "Left"
profiles["Right"] = "Right"
profiles["Fovea"] = "Fovea"
profiles["Periphery"] = "Periphery"
if "sex" in profiles:
    sex_map = {
        0: "Female",
        1: "Male"}
    profiles["sex"] = profiles["sex"].replace(sex_map)

complete_charts = []
for scalar in tqdm(["fa", "md", "mk"]):
    column_charts = []
    for column in tqdm(columns, leave=False):
        if scalar == "mk":
            x_title = f"Position ({directions_formal[column]})"
            x_labels = True
        else:
            x_title = ""
            x_labels = False
        if scalar == "fa":
            top_title = column_names_formal[column]
            if column == "LR" and "cross_fp" not in columns and not lr_only:
                top_title = top_title + ", " + column_names_formal[side_assoc[sides[0]]]
        else:
            top_title = ""
        if "profs" in dataset or column == "fov" or column == "LR_only" or column == "CST_L" or column == "UNC_L":
            y_title = scalar.upper()
            y_labels = True
        elif column == "LR" or column == "LR_fov" or (no_LR and column == "crossfp") or column == "LR_cross" or column == "LR_CST" or column == "LR_UNC":
            y_title = f"ACI, {scalar.upper()}"
            y_labels = True
        else:
            y_title = ""
            y_labels = False
        if "match" in dataset and y_title != "":
            y_title = "ACI (test+), " + y_title
        if column in ["L_OR", "R_OR", "fov", "mac", "periph", "CST", "UNC", "LR_only", "CST_L", "CST_R", "UNC_L", "UNC_R"]:
            if "match" in dataset:
                this_scale = alt.Scale(domain=aci_domain2)
            elif "cst" in b_name or "CST" in column:
                this_scale = alt.Scale(domain=cst_scalar_domains[scalar])
            elif "unc" in b_name or "UNC" in column:
                this_scale = alt.Scale(domain=unc_scalar_domains[scalar])
            elif "ukbb" in dataset or "profs" in dataset:
                this_scale = alt.Scale(domain=ukbb_scalar_domains[scalar])
            else:
                this_scale = alt.Scale(domain=scalar_domains[scalar])
        elif "match" in dataset:
            this_scale = alt.Scale(domain=match_aci_domain)
        else:
            this_scale = alt.Scale(domain=aci_domain)
        side_charts = []
        for side in sides:
            if lr_only:
                if column == "column":
                    color_name = "ACI"
                elif side == "L":
                    color_name = "Left"
                else:
                    color_name = "Right"
            elif side == "L":
                if column == "LR":
                    color_name = "Fovea"
                else:
                    color_name = "Left"
            elif side == "R":
                if column == "LR":
                    color_name = "Periphery"
                else:
                    color_name = "Right"
            else:
                if column == "LR_fov":
                    color_name = "Fovea"
                else:
                    color_name = "Periphery"
                
            this_x = alt.X('nodeID:Q', title=x_title, axis=alt.Axis(labels=x_labels))
            if f"dki_{scalar}_{get_column[column][side]}_mean" not in profiles.columns:
                raise ValueError(f"dki_{scalar}_{get_column[column][side]}_mean not in profiles.columns")
            if f"dki_{scalar}_{get_column[column][side]}_low_CI" not in profiles.columns:
                raise ValueError(f"dki_{scalar}_{get_column[column][side]}_low_CI not in profiles.columns")
            if f"dki_{scalar}_{get_column[column][side]}_high_CI" not in profiles.columns:
                raise ValueError(f"dki_{scalar}_{get_column[column][side]}_high_CI not in profiles.columns")
            this_y = alt.Y(
                f"dki_{scalar}_{get_column[column][side]}_mean:Q",
                scale=this_scale,
                title=y_title,
                axis=alt.Axis(labels=y_labels)
            )
            this_y_low_ci = alt.Y(
                f"dki_{scalar}_{get_column[column][side]}_low_CI:Q",
                scale=this_scale,
                title=y_title,
                axis=alt.Axis(labels=y_labels)
            )
            this_y_high_ci = alt.Y(
                f"dki_{scalar}_{get_column[column][side]}_high_CI:Q",
                scale=this_scale,
                title=y_title,
                axis=alt.Axis(labels=y_labels)
            )
            this_y_low_iqr = alt.Y(
                f"dki_{scalar}_{get_column[column][side]}_low_IQR:Q",
                scale=this_scale,
                title=y_title,
                axis=alt.Axis(labels=y_labels)
            )
            this_y_high_iqr = alt.Y(
                f"dki_{scalar}_{get_column[column][side]}_high_IQR:Q",
                scale=this_scale,
                title=y_title,
                axis=alt.Axis(labels=y_labels)
            )
            kwargs = {}
            uses_color = False
            if not all_sides:
                if binning is not None:
                    bin_title = binning_to_title[binning]
                    if bin_title == "Glaucoma Status":
                        if 0 in profiles[binning]:
                            profiles[binning] = profiles[binning].replace({0: "Control", 1: "Glaucoma"})
                    if not binning_to_color[binning]:
                        kwargs["strokeDash"] = alt.StrokeDash(
                            f"{binning}:N",
                            legend=alt.Legend(title=bin_title),
                            sort=binning_to_order[binning])
                    else:
                        kwargs["color"] = alt.Color(
                            f"{binning}:N",
                            legend=alt.Legend(title=bin_title),
                            #scale=alt.Scale(scheme="plasma"),
                            sort=binning_to_order[binning])
                        uses_color = True
                if not uses_color:
                    kwargs["color"] = alt.Color(
                        "age_bin",
                        legend=alt.Legend(title="Age Bin"),
                        scale=alt.Scale(scheme="plasma"),
                        sort=age_bin_order)
            elif b_name == "":
                if color_name not in profiles.columns:
                    raise ValueError(f"{color_name} not in profiles.columns")
                if "match" in dataset:
                    legend_name = "Hemisphere"
                else:
                    legend_name = "Sub-Bundle"
                kwargs["color"] = alt.Color(
                    color_name,
                    legend=alt.Legend(title=legend_name),
                    scale=alt.Scale(scheme="tableau10"),
                    sort=order)
            side_charts.append(
                alt.Chart(title=top_title).mark_line(size=line_size).encode(
                    x=this_x,
                    y=this_y,
                    **kwargs))
            if all_sides or binning is None or binning_to_color[binning]:
                if use_iqr:
                    side_charts.append(
                        alt.Chart(title=top_title).mark_line(opacity=0.5, strokeDash=[1,1]).encode(
                            x=this_x,
                            y=this_y_low_iqr,
                            **kwargs))
                    side_charts.append(
                        alt.Chart(title=top_title).mark_line(opacity=0.5, strokeDash=[1,1]).encode(
                            x=this_x,
                            y=this_y_high_iqr,
                            **kwargs))
                    side_charts.append(
                        alt.Chart(title=top_title).mark_line(opacity=0.5).encode(
                            x=this_x,
                            y=this_y_low_ci,
                            **kwargs))
                    side_charts.append(
                        alt.Chart(title=top_title).mark_line(opacity=0.5).encode(
                            x=this_x,
                            y=this_y_high_ci,
                            **kwargs))
                else:
                    side_charts.append(
                        alt.Chart(title=top_title).mark_line(opacity=0.5).encode(
                            x=this_x,
                            y=this_y_low_ci,
                            **kwargs))
                    side_charts.append(
                        alt.Chart(title=top_title).mark_line(opacity=0.5).encode(
                            x=this_x,
                            y=this_y_high_ci,
                            **kwargs))
        if ((column == "crossfp" or column == "crossfm" or "LR" in column) and column != "LR_only") or dataset == "match":
            red_line_chart = alt.Chart(title=top_title).mark_rule(color='red').encode(y='0_LINE')
            side_charts.append(red_line_chart)
        column_charts.append(alt.LayerChart(
            layer=side_charts))
    complete_charts.append(alt.HConcatChart(hconcat=column_charts))
# if uses_color:
#     profiles = profiles[profiles.age_bin == "64-67"]
profile_charts = alt.VConcatChart(vconcat=complete_charts, data=profiles)

if not all_sides:
    f_suf = f"_{sides[0]}{b_name}"
elif b_name != "":
    f_suf = f"{b_name}"
else:
    f_suf = ""

if binning is not None and binning != "is_dis":
    f_suf = f_suf + f"_{binning}"

if filt:
    f_suf = f_suf + "_filt"
if oab:
    f_suf = f_suf + "_oab"

profile_charts.configure_axis(
        labelFontSize=font_size,
        titleFontSize=font_size,
        labelLimit=0).configure_legend(
            labelFontSize=font_size,
            titleFontSize=font_size,
            titleLimit=0,
            labelLimit=0,
            columns=column_count[dataset][sides[0]],
            symbolStrokeWidth=line_size*10,
            symbolSize=line_size*100,
            orient='bottom'
        ).configure_title(
            fontSize=font_size
            ).save(f'output/first_paper_plots/tract_profiles_{f_pre}{dataset}{f_suf}.html')
