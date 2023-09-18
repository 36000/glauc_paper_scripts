import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import altair as alt
import numpy as np
import pandas as pd
import os.path as op
import os
from math import ceil, floor
from pymatch.Matcher import Matcher 
from tqdm import tqdm
import scipy.stats as st
from scipy.optimize import linear_sum_assignment
from sklearn.linear_model import LogisticRegression
from numpy.random import default_rng
import altair as alt
from sklearn.metrics import roc_auc_score, roc_curve

import compare_auc_delong_xu

font_size = 20
line_size = 3

formal_or_names = {
    "FOV": "fOR", "MAC": "mOR", "PERIPH": "pOR", "CST": "CST", "UNC": "UNC", "OR": "OR"}
b_comparisons = [
    # ("fov_L", "periph_L"),
    # ("fov_R", "periph_R"),
    # ("fov_L", "mac_L"),
    # ("fov_R", "mac_R"),
    # ("fov_L", "fov_R"),
    # ("mac_L", "mac_R"),
    # ("periph_L", "periph_R"),
    ("L_OR", "R_OR"),
    ("CST_L", "CST_R"),
    ("UNC_L", "UNC_R")]
bundle_names = [
    # "fov_L", "mac_L", "periph_L", "fov_R", "mac_R", "periph_R",
    "L_OR", "R_OR",
    "CST_L", "CST_R", "UNC_L", "UNC_R"]

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    """Modified from implementation in: """
    """https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python/33532498#33532498"""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()

    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def print_delong_table(target, hats):
    if -1 in target:
        target[target == -1] = 0
        for ii in range(len(hats)):
            hats[ii] = (hats[ii] + 1)/2

    # table_str = "      | FOV    | MAC    | PERIPH | CST    | UNC    | \n"
    # bundle_dict = {
    #     0: "FOV   ",
    #     1: "MAC   ",
    #     2: "PERIPH",
    #     3: "CST   ",
    #     4: "UNC   ",
    # }
    # comparisons = {
    #     0: {0: False, 1: True, 2: True, 3: True, 4: True},
    #     1: {0: False, 1: False, 2: True, 3: True, 4: True},
    #     2: {0: False, 1: False, 2: False, 3: True, 4: True},
    #     3: {0: False, 1: False, 2: False, 3: False, 4: False},
    #     4: {0: False, 1: False, 2: False, 3: False, 4: False},
    # }
    table_str = "      | OR     | CST    | UNC    | \n"
    bundle_dict = {
        0: "OR    ",
        1: "CST   ",
        2: "UNC   ",
    }
    comparisons = {
        0: {0: False, 1: True, 2: True},
        1: {0: False, 1: False, 2: False},
        2: {0: False, 1: False, 2: False},
    }

    pvals = []
    for ii in range(len(hats)):
        for jj in range(len(hats)):
            if comparisons[ii][jj]:
                pvals.append(10**compare_auc_delong_xu.delong_roc_test(
                    target, hats[ii], hats[jj])[0][0])

    import statsmodels
    _, pvals = statsmodels.stats.multitest.fdrcorrection(pvals)

    pval_index = 0
    for ii in range(3):
        bundle = bundle_dict[ii]
        table_str += f"{bundle}| "
        for jj in range(len(hats)):
            if comparisons[ii][jj]:
                pval = pvals[pval_index]
                pval_index = pval_index + 1
                if np.isnan(pval):
                    table_str += f"N/A    | "
                else:
                    table_str += f"{'%.4f' % round(pval, 4)} | "
            else:
                table_str += f"N/A    | "
        table_str += "\n"
    print(table_str)


    
    # table_str = "FOV   | MAC   | PERIPH| CST   | UNC   | \n"
    # for ii in range(len(hats)):
    #     pval = st.mannwhitneyu(hats[ii], target, use_continuity=False, alternative="greater").pvalue
    #     print(len(hats[ii]))
    #     table_str += f"{'%.3f' % round(pval, 3)} | "
    # print(table_str)
        

def auc_score(target, hat):
    auc, ci = get_auc_unc(target, hat)
    return auc-ci, auc, auc+ci

def get_auc_unc(target, hat):
    if -1 in target:
        target[target == -1] = 0
        hat = (hat + 1)/2
    auc, var = compare_auc_delong_xu.delong_roc_variance(target, hat)
    dunn_correct = 5
    ci_fac = st.norm.ppf(1-(0.05/dunn_correct)/2)
    return auc, ci_fac * np.sqrt(var)


def make_roc_plot(fname, roc_df, auc_df, b_order, b_order_auc, title):
    baseline = alt.Chart(roc_df).mark_line(size=line_size, strokeDash=[20,5], color = 'black').encode(
        alt.X('dotted_line', title='FPR'),
        alt.Y('dotted_line', title='TPR'))
    baseline.properties(title=title).configure_axis(
        labelFontSize=font_size,
        titleFontSize=font_size
        ).configure_legend(
            labelFontSize=font_size,
            titleFontSize=font_size,
            labelLimit=0,
            symbolStrokeWidth=10,
        ).configure_title(
            fontSize=font_size
            ).save(f"output/pheno_match/{fname}_baseline.html")

    roc_df_ctrl = roc_df[roc_df["Bundle"].isin(b_order[-2:])]
    roc_line = alt.Chart(roc_df_ctrl).mark_line(size=line_size).encode(
        x='FPR',
        y='TPR',
        color=alt.Color('Bundle',
            sort=b_order, scale=alt.Scale(
                domain=b_order,
                #range=['#750404','#D12828', '#FA6161', '#071bf5', '#626efb']))
                range=['#750404', '#071bf5', '#626efb']))
    )

    chart = (roc_line + baseline).properties(
        title=title)
    chart.configure_axis(
        labelFontSize=font_size,
        titleFontSize=font_size
        ).configure_legend(
            labelFontSize=font_size,
            titleFontSize=font_size,
            labelLimit=0,
            symbolStrokeWidth=10,
        ).configure_title(
            fontSize=font_size
            ).save(f"output/pheno_match/{fname}_ctrl.html")

    roc_line = alt.Chart(roc_df).mark_line(size=line_size).encode(
        x='FPR',
        y='TPR',
        color=alt.Color('Bundle',
            sort=b_order, scale=alt.Scale(
                domain=b_order,
                #range=['#750404','#D12828', '#FA6161', '#071bf5', '#626efb']))
                range=['#750404', '#071bf5', '#626efb']))
    )

    chart = (roc_line + baseline).properties(
        title=title)
    base = alt.Chart(auc_df).encode(
        x=alt.X("Bundle", sort=b_order_auc)
    )
    rule = base.mark_rule(color = 'black', size=line_size).encode(
        y=alt.Y("lci", title=None),
        y2=alt.Y2("hci", title=None)
    )
    auc_chart = base.mark_bar().encode(
        y=alt.Y("AUC", scale=alt.Scale(domain=[0.4, 0.8])),
        color=alt.Color(
            'Bundle',
            sort=b_order_auc, scale=alt.Scale(
                domain=b_order_auc,
                range=['#750404','#D12828', '#FA6161', '#071bf5', '#626efb']),
            legend=alt.Legend(values=b_order)),
    )
    # auc_chart = auc_chart + rule
    # chart = chart | auc_chart
    chart.configure_axis(
        labelFontSize=font_size,
        titleFontSize=font_size
        ).configure_legend(
            labelFontSize=font_size-5,
            titleFontSize=font_size-5,
            labelLimit=0,
            symbolStrokeWidth=10,

            columns=2,
            orient='right'
        ).configure_title(
            fontSize=font_size
            ).save(f"output/pheno_match/{fname}.html")


def match_pipe():
    pheno = pd.read_csv("output/pheno.csv", low_memory=False)
    if "glauc" in pheno_of_interest:
        dis_subs = pheno[np.logical_not(
            np.isnan(pheno["6119-0.0"]) &
            np.isnan(pheno["6119-1.0"]) &
            np.isnan(pheno["6119-2.0"]) &
            np.isnan(pheno["6119-3.0"]))].eid.unique()
    elif pheno_of_interest == "logmar":
        dis_subs = pheno[(
            ((pheno["5208-0.0"] <= 0.3) & (pheno["5201-0.0"] <= 0.3)) |
            ((pheno["5208-1.0"] <= 0.3) & (pheno["5201-1.0"] <= 0.3)))].eid.unique()
    elif pheno_of_interest == "ambly":
        dis_subs = pheno[np.logical_not(
            np.isnan(pheno["5408-0.0"]) &
            np.isnan(pheno["5408-1.0"]) &
            np.isnan(pheno["5408-2.0"]) &
            np.isnan(pheno["5408-3.0"]))].eid.unique()
    elif pheno_of_interest == "any":
        dis_subs = pheno[
            (pheno["6148-2.0"] != -7)
            & (pheno["6148-3.0"] != -7)].eid.unique()
    # elif pheno_of_interest == "alzh":
    #     dis_subs = pheno[].eid.unique()
    elif "age" in pheno_of_interest:
        for_dis_pheno = pheno
        if "nocut" not in pheno_of_interest:
            for acuity_f in ["5201-0.0", "5201-1.0", "5208-0.0", "5208-1.0"]:
                for_dis_pheno = for_dis_pheno[(for_dis_pheno[acuity_f] <= 0.3) | np.isnan(for_dis_pheno[acuity_f])]
        for_dis_pheno = for_dis_pheno[
            (for_dis_pheno["6148-2.0"] == -7)
            | (for_dis_pheno["6148-3.0"] == -7)]
        for_dis_pheno = for_dis_pheno[
                np.isnan(for_dis_pheno["6119-0.0"]) &
                np.isnan(for_dis_pheno["6119-1.0"]) &
                np.isnan(for_dis_pheno["6119-2.0"]) &
                np.isnan(for_dis_pheno["6119-3.0"])]
        dis_subs = for_dis_pheno[(for_dis_pheno["21003-2.0"] >= 70) & (for_dis_pheno["21003-2.0"] < 80)].eid.unique()
    else:
        raise ValueError("unknown pheno_of_interest")

    #print(len(pheno))
    # exclude acuity worse than 0.3
    if "nocut" not in pheno_of_interest:
        for acuity_f in ["5201-0.0", "5201-1.0", "5208-0.0", "5208-1.0"]:
            pheno = pheno[(pheno[acuity_f] <= 0.3) | np.isnan(pheno[acuity_f])]

    # exclude glauc
    pheno = pheno[
        (np.isnan(pheno["6119-0.0"]) &
         np.isnan(pheno["6119-1.0"]) &
         np.isnan(pheno["6119-2.0"]) &
         np.isnan(pheno["6119-3.0"])) | pheno.eid.isin(dis_subs)]

    # exclude other retinal disease
    pheno = pheno[
        (pheno["6148-2.0"] == -7)
        | (pheno["6148-3.0"] == -7)
        | pheno.eid.isin(dis_subs)]
    #print(len(pheno))

    pheno = pheno.rename(columns=pheno_dict)
    pheno_all = pheno[["eid", *all_pheno_names]]
    pheno = pheno[["eid", *pheno_names]]

    test = pheno[pheno.eid.isin(dis_subs)]
    control = pheno[~pheno.eid.isin(dis_subs)]
    if "age" in pheno_of_interest:
        control["age"] = control["age"] + 10
    test_all = pheno_all[pheno_all.eid.isin(dis_subs)]
    control_all = pheno_all[~pheno_all.eid.isin(dis_subs)]
    if "age" in pheno_of_interest:
        control_all["age"] = control_all["age"] + 10

    test.loc[:, f'{pheno_of_interest}_status'] = 1
    control.loc[:, f'{pheno_of_interest}_status'] = 0
    test_all.loc[:, f'{pheno_of_interest}_status'] = 1
    control_all.loc[:, f'{pheno_of_interest}_status'] = 0
    print(test)
    print(control)
    print(f"Test mean+/-sd age before match: {np.nanmean(test_all['age'])}+/-{np.nanstd(test_all['age'])}")
    print(f"Control mean+/-sd age before match: {np.nanmean(control_all['age'])}+/-{np.nanstd(control_all['age'])}")

    if "none" in pheno_of_interest:
        test_eid = test.eid.to_numpy()
        ctrl_eid = np.random.choice(control.eid.to_numpy(), len(test_eid), replace=False)

        np.savetxt(
            f"output/pheno_match/{pheno_of_interest}/ctrl_sub.txt",
            ctrl_eid)
        match_df = pd.DataFrame()
        for eids, frame in zip([test_eid, ctrl_eid], [test_all, control_all]):
            for matchid, eid in enumerate(eids):
                _dict = frame[frame.eid==eid].iloc[0].to_dict()
                _dict["match_id"] = matchid
                match_df = match_df.append(_dict, ignore_index=True)
        for intcolumn in ["match_id", "eid", f'{pheno_of_interest}_status']:
            match_df[intcolumn] = match_df[intcolumn].astype(int)
        print("Len of matched eid: " + str(len(test_eid)))
        match_df.to_csv(f"output/pheno_match/{pheno_of_interest}/match.csv", index=False)
    elif matcher == "psm":
        # for reproducibility
        np.random.seed(20170925)
        test = test.fillna(-1)
        control = control.fillna(-1)
        m = Matcher(test, control, yvar=f'{pheno_of_interest}_status', exclude=["eid"])
        m.fit_scores(balance=True, nmodels=100)
        m.match(method="min", nmatches=1, threshold=0.001, with_replacement=False)
        m.record_frequency()
        print(m.record_frequency())
        print(m.matched_data)
        print("Len of matched eid: " + str(len(m.matched_data.eid)))
        print("Len of matched unique eid: " + str(len(m.matched_data.eid.unique())))

        m.matched_data.replace(-1, np.nan).to_csv(f"output/pheno_match/{pheno_of_interest}/match.csv", index=False)
        np.savetxt(
            f"output/pheno_match/{pheno_of_interest}/ctrl_sub.txt",
            m.matched_data[m.matched_data[f'{pheno_of_interest}_status'] == 0].eid)
    elif matcher == "mdm":
        test_wo_eid = test.drop(["eid", f'{pheno_of_interest}_status'], axis=1).to_numpy()
        print("Original test len: " + str(test_wo_eid.shape[0]))
        ctrl_wo_eid = control.drop(["eid", f'{pheno_of_interest}_status'], axis=1).to_numpy()
        print("Original ctrl len: " + str(ctrl_wo_eid.shape[0]))
        test_wo_eid[np.isnan(test_wo_eid)] = -1  # TODO: try this as mean
        ctrl_wo_eid[np.isnan(ctrl_wo_eid)] = -1  # TODO: correlate sibjects with missing scans with other interesting variables

        v_inv = np.linalg.inv(np.cov(np.concatenate((test_wo_eid, ctrl_wo_eid), axis=0).T, ddof=0))
        diff = test_wo_eid[:, np.newaxis] - ctrl_wo_eid[np.newaxis, :]
        nbrs = np.sqrt(np.sum((diff@v_inv)*diff, axis=2))

        # nbrs2 = np.zeros((test_wo_eid.shape[0], ctrl_wo_eid.shape[0]))
        # for i in range(test_wo_eid.shape[0]):
        #     for j in range(ctrl_wo_eid.shape[0]):
        #         nbrs2[i, j] = mahalanobis(test_wo_eid[i], ctrl_wo_eid[j], v_inv)

        row_ind, col_ind = linear_sum_assignment(nbrs)
        threshold = 0.3
        if pheno_of_interest == "glauc_sec":
            threshold = 0.3
        if pheno_of_interest == "glauc_mega":
            threshold = 0.5
        filtered_row_ind = [] # remove matches that are too bad, 0.2
        filtered_col_ind = []
        for row_idx, col_idx in zip(row_ind, col_ind):
            if nbrs[row_idx, col_idx] < threshold:
                filtered_row_ind.append(row_idx)
                filtered_col_ind.append(col_idx)
        test_eid = test.eid.to_numpy()[filtered_row_ind]
        ctrl_eid = control.eid.to_numpy()[filtered_col_ind]

        np.savetxt(
            f"output/pheno_match/{pheno_of_interest}/ctrl_sub.txt",
            ctrl_eid)
        match_df = pd.DataFrame()
        for eids, frame in zip([test_eid, ctrl_eid], [test_all, control_all]):
            for matchid, eid in enumerate(eids):
                _dict = frame[frame.eid==eid].iloc[0].to_dict()
                _dict["match_id"] = matchid
                match_df = match_df.append(_dict, ignore_index=True)
        for intcolumn in ["match_id", "eid", f'{pheno_of_interest}_status']:
            match_df[intcolumn] = match_df[intcolumn].astype(int)
        print("Len of matched eid: " + str(len(test_eid)))

        test_post_match = match_df[match_df[f'{pheno_of_interest}_status'].astype(bool)]
        ctrl_post_match = match_df[~match_df[f'{pheno_of_interest}_status'].astype(bool)]
        print(f"Test mean+/-sd age after match: {np.nanmean(test_post_match['age'])}+/-{np.nanstd(test_post_match['age'])}")
        print(f"Control mean+/-sd age after match: {np.nanmean(ctrl_post_match['age'])}+/-{np.nanstd(ctrl_post_match['age'])}")
        match_df.to_csv(f"output/pheno_match/{pheno_of_interest}/match.csv", index=False)


def wide_pipe():
    dis_df = pd.read_csv(f"output/pheno_match/{pheno_of_interest}/match.csv")
    dis_df = dis_df.drop(columns=["scores", "record_id"], errors="ignore")

    dis_df_wide = dis_df.pivot(
        index="match_id", columns=f'{pheno_of_interest}_status',
        values=[*all_pheno_names, "eid"]).reset_index()
    dis_df_wide.columns = dis_df_wide.columns.map(str).map(''.join).astype(str).str.replace("'", "")
    print(dis_df_wide)
    dis_df_wide = dis_df_wide.drop_duplicates(["(eid, 0)"])
    print(dis_df_wide)
    dis_df_wide.to_csv(f"output/pheno_match/{pheno_of_interest}/wide.csv", index=False)


def ase_info_pipe():
    dis_df_wide = pd.read_csv(f"output/pheno_match/{pheno_of_interest}/wide.csv")
    charts = []
    for match_val in all_pheno_names:
        one = dis_df_wide[f"({match_val}, 0)"].to_numpy()
        one[np.isnan(one)] = -1
        two = dis_df_wide[f"({match_val}, 1)"].to_numpy()
        two[np.isnan(two)] = -1
        _min = min(one[one != -1].min(), two[two != -1].min())
        _max = max(one.max(), two.max())
        if match_val == "age":
            _bins = np.arange(_min, _max)
        elif match_val == "sex":
            _bins = 2
        elif match_val == "ethnicity":
            _bins = np.unique(np.concatenate((one, two)))
        elif match_val == "FVC":
            _bins = np.concatenate(([-1.5, -0.5], np.arange(floor(_min), ceil(_max), step=0.5)))
        elif match_val == "PEF":
            _bins = np.concatenate(([-1.5, -0.5], np.arange(floor(_min/50)*50, ceil(_max/50)*50, step=50)))
        elif match_val == "TDI":
            _bins = np.arange(-7, 11, step=1)
        else:
            raise ValueError(f"bin algo not found for {match_val}")
        H, xedges, yedges = np.histogram2d(
            one, two,
            bins=_bins, density=False)
        x, y = np.meshgrid(xedges[:-1], yedges[:-1])
        this_data = pd.DataFrame({
            'x': x.ravel(),
            'y': y.ravel(),
            match_val: H.ravel()})
        charts.append(alt.Chart(this_data).mark_rect().encode(
            x=f"x:O",
            y=f"y:O",
            color=alt.Color(match_val, scale=alt.Scale(type="symlog")),
            tooltip=[match_val]
        ).interactive())

    scatters = alt.HConcatChart(hconcat=charts)


    charts = []
    for match_val in all_pheno_names:
        dis_df_wide[f"{match_val}_matched"] = (dis_df_wide[f"({match_val}, 0)"] == dis_df_wide[f"({match_val}, 1)"])
        charts.append(alt.Chart(dis_df_wide).mark_bar().encode(
            x=f"{match_val}_matched:Q",
            y=alt.Y(f"count({match_val}_matched):Q", scale=alt.Scale(type="symlog")),
            tooltip=[f"count({match_val}_matched)"]
        ).interactive())

    counts = alt.HConcatChart(hconcat=charts)

    charts = []
    steps = {
        "age": 3,
        "sex": 0.5,
        "ethnicity": 500,
        "FVC": 0.5,
        "PEF": 50,
        "TDI": 1
    }
    for match_val in all_pheno_names:
        to_layer = []
        for ii, color in zip(range(2), ["blue", "red"]):  # 0 is control
            to_layer.append(alt.Chart(dis_df_wide, title=match_val).mark_bar(color=color, opacity=0.6).encode(
                alt.X(
                    f"({match_val}, {ii})",
                    bin=alt.Bin(step=steps[match_val]),
                    title=match_val),
                y='count()'))
        charts.append(alt.LayerChart(layer=to_layer).interactive())

    hist = alt.HConcatChart(hconcat=charts)

    (scatters & counts & hist).configure_axis(
        labelFontSize=font_size,
        titleFontSize=font_size
        ).configure_legend(
            labelFontSize=font_size,
            titleFontSize=font_size,
            labelLimit=0,
        ).configure_title(
            fontSize=font_size
            ).save(f"output/pheno_match/{pheno_of_interest}/ase_info.html")

    print("perfect matches: " + str(np.sum(
        (dis_df_wide[f"(age, 0)"] == dis_df_wide[f"(age, 1)"])
        & (dis_df_wide[f"(sex, 0)"] == dis_df_wide[f"(sex, 1)"])
        & (dis_df_wide[f"(ethnicity, 0)"] == dis_df_wide[f"(ethnicity, 1)"])
    )))
    print("out of: " + str(len(dis_df_wide["(age, 0)"])))


def profiles_pipe():
    profiles = pd.read_csv("output/tract_profiles_wide.csv")
    print(len(dis_df_wide["(eid, 0)"].unique()))
    print(len(dis_df_wide["(eid, 1)"].unique()))
    matched_profiles = pd.DataFrame()
    for _, row in tqdm(dis_df_wide.iterrows()):
        ctrl_profile = profiles[profiles.subjectID == row["(eid, 0)"]].reset_index(drop=True)
        dis_profile = profiles[profiles.subjectID == row["(eid, 1)"]].reset_index(drop=True)
        comb_profile = ctrl_profile.join(
            dis_profile, on="nodeID", how='outer',
            lsuffix="_ctrl", rsuffix="_dis", sort=False)
        comb_profile = comb_profile.assign(match_id=row["(match_id, )"])
        matched_profiles = matched_profiles.append(
            comb_profile, ignore_index=True)
    #print(matched_profiles)
    print(len(matched_profiles.subjectID_ctrl.unique()))
    #print(len(matched_profiles))
    print(len(matched_profiles))

    matched_profiles.to_csv(f"output/pheno_match/{pheno_of_interest}/profiles.csv")


def aci_pipe():
    dataframe = pd.read_csv(f"output/pheno_match/{pheno_of_interest}/profiles.csv")
    dataframe = dataframe[(dataframe.nodeID >= 20) & (dataframe.nodeID < 80)]

    new_bundles = []
    for scalar in tqdm(["dki_fa", "dki_md", "dki_mk"]):
        for bundle in bundle_names:
            column_a = dataframe[f"{scalar}_{bundle}_dis"]
            column_b = dataframe[f"{scalar}_{bundle}_ctrl"]
            dataframe[f"{scalar}_{bundle}"] = 2*(
                column_a-column_b)/(np.abs(column_a)+np.abs(column_b))
        for bundle_a, bundle_b in b_comparisons:
            column_a = dataframe[f"{scalar}_{bundle_a}"]
            column_b = dataframe[f"{scalar}_{bundle_b}"]
            dataframe[f"{scalar}_{bundle_a}_{bundle_b}_aci"] = 2*(
                column_a-column_b)/(np.abs(column_a)+np.abs(column_b))
            if scalar == "dki_fa":
                new_bundles.append(f"{bundle_a}_{bundle_b}_aci")
    bundles = [*bundle_names, *new_bundles]

    def get_mean_ci(this_bundle):
        this_mean = np.nanmean(this_bundle)
        this_low, this_high = st.bootstrap(
            (this_bundle,), np.nanmean, axis=0,
            confidence_level=0.95, n_resamples=10000,
            method="BCa").confidence_interval
        this_low_iqr = np.nanpercentile(this_bundle, 25)
        this_high_iqr = np.nanpercentile(this_bundle, 75)
        return this_mean, this_low, this_high, this_low_iqr, this_high_iqr

    new_dataframe = pd.DataFrame()
    dataframe_profs = pd.DataFrame()
    for nodeID in tqdm(range(20, 80)):
        this_row = {"nodeID": nodeID}
        this_row_profs = {"nodeID": nodeID, "is_dis": 1}
        this_row_profs_ctrl = {"nodeID": nodeID, "is_dis": 0}
        this_bundle = dataframe[dataframe.nodeID == nodeID]
        for column in tqdm(bundle_names, leave=False):
            for scalar in tqdm(["dki_fa", "dki_md", "dki_mk"], leave=False):
                c_name = f"{scalar}_{column}"         
                this_row[f"{c_name}_mean"], this_row[f"{c_name}_low_CI"], this_row[f"{c_name}_high_CI"], this_row[f"{c_name}_low_IQR"], this_row[f"{c_name}_high_IQR"] =\
                    get_mean_ci(this_bundle[c_name])
                if "aci" not in column:
                    this_row_profs[f"{c_name}_mean"], this_row_profs[f"{c_name}_low_CI"], this_row_profs[f"{c_name}_high_CI"], this_row_profs[f"{c_name}_low_IQR"], this_row_profs[f"{c_name}_high_IQR"] =\
                        get_mean_ci(this_bundle[c_name + "_dis"])
                    this_row_profs_ctrl[f"{c_name}_mean"], this_row_profs_ctrl[f"{c_name}_low_CI"], this_row_profs_ctrl[f"{c_name}_high_CI"], this_row_profs_ctrl[f"{c_name}_low_IQR"], this_row_profs_ctrl[f"{c_name}_high_IQR"] =\
                        get_mean_ci(this_bundle[c_name + "_ctrl"])
        new_dataframe = new_dataframe.append(this_row, ignore_index=True)
        dataframe_profs = dataframe_profs.append(this_row_profs, ignore_index=True)
        dataframe_profs = dataframe_profs.append(this_row_profs_ctrl, ignore_index=True)

    new_dataframe.to_csv(f"output/aci_for_paper/profiles_w_aci_{pheno_of_interest}_match.csv", index=False)
    dataframe_profs.to_csv(f"output/aci_for_paper/profiles_w_aci_{pheno_of_interest}_profs.csv", index=False)


def fracridge_pipe():
    os.makedirs(f"output/pheno_match/{pheno_of_interest}/fracridge_derivs/", exist_ok=True)
    dataframe = pd.read_csv(f"output/pheno_match/{pheno_of_interest}/profiles.csv")

    all_subjects = dataframe.subjectID_ctrl.unique()
    print(len(all_subjects))

    rng = default_rng(12345)
    indices = rng.permutation(len(all_subjects))
    cut1 = int(0.8*len(all_subjects))
    training_idx, test_idx = indices[:cut1], indices[cut1:]
    cut2 = int(0.8*len(training_idx))
    training_idx, val_idx = indices[:cut2], indices[cut2:cut1]
    subs = {
        "train": all_subjects[training_idx],
        "test": all_subjects[test_idx],
        "val": all_subjects[val_idx]}
    # bd = {
    #     "FOV": ["fov_L",  "fov_R"],
    #     "MAC": ["mac_L",  "mac_R"],
    #     "PERIPH": ["periph_L", "periph_R"],
    #     "CST": ["CST_L", "CST_R"],
    #     "UNC": ["UNC_L", "UNC_R"]}
    bd = {
        "OR": ["L_OR",  "R_OR"],
        "CST": ["CST_L", "CST_R"],
        "UNC": ["UNC_L", "UNC_R"]}
    Xs = {}
    already_exists = True
    for b_name, bundles in bd.items():
        Xs[b_name] = {}
        for sett in subs.keys():
            if not op.exists((
                    f"output/pheno_match/{pheno_of_interest}/"
                    f"fracridge_derivs/X_{b_name}_{sett}.npy")):
                already_exists = False
            Xs[b_name][sett] = np.zeros((len(subs[sett]), 2, len(bundles), 3, 100))
    if False or global_recompute or not already_exists:
        for sett in tqdm(subs.keys(), leave=False):
            for ii, subjectID in enumerate(tqdm(subs[sett], leave=False)):
                sub_df = dataframe[dataframe.subjectID_ctrl == subjectID]
                for nodeID in tqdm(range(100), leave=False):
                    node_df = sub_df[sub_df.nodeID == nodeID]
                    for b_name, bundles in tqdm(bd.items(), leave=False):
                        for ll, bundle in enumerate(tqdm(bundles, leave=False)):
                            for jj, scalar in enumerate(tqdm(["dki_fa", "dki_md", "dki_mk"], leave=False)):
                                for kk, suf in enumerate(tqdm(["ctrl", "dis"], leave=False)):
                                    Xs[b_name][sett][ii, kk, ll, jj, nodeID] = node_df[
                                        f"{scalar}_{bundle}_{suf}"].to_numpy()[0]
        for b_name in bd.keys():
            for sett in subs.keys():
                np.save((
                    f"output/pheno_match/{pheno_of_interest}/"
                    f"fracridge_derivs/X_{b_name}_{sett}.npy"), Xs[b_name][sett])
    else:
        for b_name in bd.keys():
            for sett in subs.keys():
                Xs[b_name][sett] = np.load((
                    f"output/pheno_match/{pheno_of_interest}/"
                    f"fracridge_derivs/X_{b_name}_{sett}.npy"))

    # remove nans
    # c_valid = 0
    for b_name, bundles in bd.items():
        for sett in subs.keys():
            # for subject in range(Xs[b_name][sett].shape[0]):
            #     for ll, bundle in enumerate(bundles):
            #         for ii in range(2):
            #             if np.sum(np.isnan(Xs[b_name][sett][subject, ii, ll].flatten())) == 0:
            #                 c_valid = c_valid + 1
            # raise ValueError(f"{c_valid}, {Xs[b_name][sett].shape[0]*4}")
            Xs[b_name][sett][
                np.any(np.isnan(Xs[b_name][sett]), axis=(2, 3, 4))] =\
                np.nanmean(Xs[b_name][sett], axis=(0, 1))

    step_size = 2
    roc_df = pd.DataFrame()
    auc_df = pd.DataFrame()
    b_order = []
    b_order_auc = []
    hats = []
    for b_name, bundles in tqdm(bd.items()):
        this_x = Xs[b_name]["train"][:, :, :, :, 20:80:step_size]
        train_len = this_x.shape[0]

        this_x = this_x.reshape(train_len*2, -1)
        this_std = np.std(this_x, axis=0)
        this_mean = np.mean(this_x, axis=0)
        this_x = (this_x-this_mean)/this_std
        this_y = np.zeros((train_len, 2))
        this_y[:, 1] = 1
        this_y = this_y.flatten()
        lr = LogisticRegression(solver="liblinear")
        # fr = FracRidgeRegressor(jit=False)
        lr.fit(this_x, this_y)


        this_x = Xs[b_name]["test"][:, :, :, :, 20:80:step_size]
        val_len = this_x.shape[0]
        this_x = this_x.reshape(val_len*2, -1)
        this_std = np.std(this_x, axis=0)
        this_mean = np.mean(this_x, axis=0)
        this_x = (this_x-this_mean)/this_std
        target = np.zeros((val_len, 2))
        target[:, 1] = 1
        target = target.flatten()

        target_hat = lr.predict_proba(this_x)[:, 1]

        print(f"{b_name} # of subs in val: " + str(val_len))
        print(f"{b_name} logistic AUC: " + str(auc_score(target.flatten(), target_hat.flatten())))

        # fpr, tpr = roc_curve(
        #     target.flatten(),
        #     target_hat[..., 4].flatten(),
        #     np.arange(-1, 1, 0.1))
        fpr, tpr, thresh = roc_curve(
            target.flatten(),
            target_hat.flatten())
        this_roc_df = pd.DataFrame()
        this_roc_df['FPR'] = fpr
        this_roc_df['TPR'] = tpr
        this_roc_df['dotted_line'] = np.asarray(tpr > 0.5).astype(float)
        bundle_name_w_roc = formal_or_names[b_name] + f", AUC: {round(roc_auc_score(target.flatten(), target_hat.flatten()), 2)}"
        if "UNC" in bundle_name_w_roc or "CST" in bundle_name_w_roc:
            bundle_name_w_roc = "Control: " + bundle_name_w_roc
        else:
            bundle_name_w_roc = "OR: " + bundle_name_w_roc

        b_order.append(bundle_name_w_roc)
        this_roc_df['Bundle'] = bundle_name_w_roc
        roc_df = pd.concat((roc_df, this_roc_df), ignore_index=True)

        lci, auc, hci = auc_score(target.flatten(), target_hat.flatten())
        b_order_auc.append(b_name)
        auc_df = auc_df.append({"Bundle": b_name, "lci": lci, "AUC": auc, "hci": hci}, ignore_index=True)
        hats.append(target_hat.flatten())
    
    print_delong_table(target.flatten(), hats)

    make_roc_plot(
        f"{pheno_of_interest}/fracridge_roc",
        roc_df, auc_df, b_order, b_order_auc, "Logistic Regression ROC")


def cnn_pipe():
    os.makedirs(f"output/pheno_match/{pheno_of_interest}/{use_model}_models/", exist_ok=True)
    if use_model == "resnet":
        from afqinsight.nn.tf_models import cnn_resnet
        model = cnn_resnet
        model_name = "1D Convolutional Resnet"
    elif use_model == "blstm":
        from afqinsight.nn.tf_models import blstm1
        model = blstm1
        model_name = "BLSTM 1"
    else:
        raise ValueError("use_model param not recognized")
    from afqinsight.augmentation import jitter, time_warp, scaling
    import tensorflow as tf
    from tensorflow.keras.callbacks import ModelCheckpoint
    import shap as shaply
    subs = ["train", "test", "val"]
    # bd = {
    #     "FOV": ["fov_L",  "fov_R"],
    #     "MAC": ["mac_L",  "mac_R"],
    #     "PERIPH": ["periph_L", "periph_R"],
    #     "CST": ["CST_L", "CST_R"],
    #     "UNC": ["UNC_L", "UNC_R"]}
    bd = {
        "OR": ["L_OR",  "R_OR"],
        "CST": ["CST_L", "CST_R"],
        "UNC": ["UNC_L", "UNC_R"]}
    datasets = {}
    raw_dataset = {}
    batch_size = 128

    for b_name in bd.keys():
        datasets[b_name] = {}
        raw_dataset[b_name] = {}
        for sett in subs:
            this_X = np.load((
                f"output/pheno_match/{pheno_of_interest}/"
                f"fracridge_derivs/X_{b_name}_{sett}.npy"))

            target = np.zeros((this_X.shape[0], 2))
            target[:, 1] = 1
            target = target.flatten()

            this_X[
                np.any(np.isnan(this_X), axis=(2, 3, 4))] =\
                np.nanmean(this_X, axis=(0, 1))

            this_X = this_X[:, :, :, :, 10:90]
            this_X = this_X.reshape((-1, 6, 80))
            this_X = np.swapaxes(this_X, 1, 2)

            raw_dataset[b_name][sett] = this_X

            this_std = np.std(this_X, axis=0)
            this_mean = np.mean(this_X, axis=0)
            this_X = (this_X-this_mean)/this_std

            dataset = tf.data.Dataset.from_tensor_slices((this_X, target))
            dataset = dataset.batch(batch_size)
            datasets[b_name][sett] = dataset

    results = {}
    for b_name in tqdm(bd.keys()):
        if True or global_recompute or not op.exists(f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}.h5"):
            # model from the tf_models module
            model_cnn = model(
                input_shape=(80, 6) , n_classes=1, output_activation=None, verbose=True)
            # Compile the model and fit it using training data
            mcnn_save = ModelCheckpoint(
                f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}.h5",
                save_best_only=True, monitor='mean_squared_error', mode='min')
            model_cnn.compile(
                loss= "mean_squared_error" , optimizer="adam", metrics=["mean_squared_error"])
            model_cnn.fit(datasets[b_name]["train"], callbacks=[mcnn_save], epochs=300)

    def get_X_from_dataset(dataset):
        dataset = list(dataset.unbatch().as_numpy_iterator())
        this_x = np.zeros((len(dataset), 80, 6))
        this_y = np.zeros((len(dataset)))
        for i, sub_data in enumerate(dataset):
            this_x[i, :, :] = sub_data[0]
            this_y[i] = sub_data[1]
        return this_x, this_y

    roc_df = pd.DataFrame()
    auc_df = pd.DataFrame()
    b_order = []
    b_order_auc = []
    hats = []
    for b_name in bd.keys():
        model_cnn = tf.keras.models.load_model(
            f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}.h5")
        target = np.concatenate([y for x, y in datasets[b_name]["test"]], axis=0)
        target_hat = np.asarray(model_cnn.predict(datasets[b_name]["test"])[:, 0])
        print(f"{b_name}: {auc_score(target, target_hat)}")
        fpr, tpr, tresh = roc_curve(target, target_hat)
        # fpr, tpr = roc_curve(
        #     target,
        #     target_hat,
        #     np.arange(0, 1, 0.05),
        #     use_zero=True)
        this_roc_df = pd.DataFrame()
        this_roc_df['FPR'] = fpr
        this_roc_df['TPR'] = tpr
        this_roc_df['dotted_line'] = np.asarray(tpr > 0.5).astype(float)
        bundle_name_w_roc = formal_or_names[b_name] + f", AUC: {round(roc_auc_score(target, target_hat), 2)}"
        if "UNC" in bundle_name_w_roc or "CST" in bundle_name_w_roc:
            bundle_name_w_roc = "Control: " + bundle_name_w_roc
        else:
            bundle_name_w_roc = "OR: " + bundle_name_w_roc
        b_order.append(bundle_name_w_roc)
        this_roc_df['Bundle'] = bundle_name_w_roc
        roc_df = pd.concat((roc_df, this_roc_df), ignore_index=True)

        lci, auc, hci = auc_score(target, target_hat)
        b_order_auc.append(b_name)
        auc_df = auc_df.append({"Bundle": b_name, "lci": lci, "AUC": auc, "hci": hci}, ignore_index=True)
        hats.append(target_hat)

        tf.keras.utils.plot_model(
            model_cnn,
            to_file=f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}_model.png",
            show_layer_names=False,
            show_shapes=False)

    print_delong_table(target, hats)

    make_roc_plot(
        f"{pheno_of_interest}/cnn_roc",
        roc_df, auc_df, b_order, b_order_auc, f'{model_name} ROC')

    feature_names = []
    for nodeID in range(10, 90):
        for bunsca in ["FA_L", "MD_L", "MK_L", "FA_R", "MD_R", "MK_R"]:
            feature_names.append(bunsca + "_" + str(nodeID))
    feature_names = np.asarray(feature_names)
    grouped_feature_names = []
    for section in ["ANT", "POST"]:
        for bunsca in ["FA_L", "MD_L", "MK_L", "FA_R", "MD_R", "MK_R"]:
            grouped_feature_names.append(bunsca + "_" + section)
    grouped_feature_names = np.asarray(grouped_feature_names)

    for b_name, bundles in tqdm(bd.items()):
        if True or global_recompute or not op.exists(f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}_shap.npy"):
            model_cnn = tf.keras.models.load_model(
                f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}.h5")
            train_x, _ = get_X_from_dataset(datasets[b_name]["train"])
            val_x, val_y = get_X_from_dataset(datasets[b_name]["test"])
            expl = shaply.GradientExplainer(model_cnn, train_x)
            expected_value = tf.reduce_mean(expl.explainer.model(expl.explainer.data), 0)
            shap_values = expl.shap_values(val_x)

            shap_values = np.asarray(shap_values[0])
            np.save(
                f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}_shap.npy",
                shap_values)
            np.save(
                f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}_ev.npy",
                expected_value)

    for b_name, bundles in tqdm(bd.items()):
        if True:
            val_x, val_y = get_X_from_dataset(datasets[b_name]["test"])
            shap_values = np.load(f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}_shap.npy")
            val_y = val_y.astype(bool)
            shap_glauc = shap_values[val_y]
            shap_ctrl = shap_values[~val_y]
            val_x_glauc = val_x[val_y]
            val_x_ctrl = val_x[~val_y]

            shaps = shap_values.reshape((-1, 80*6))

            shaps_grouped = shap_values.reshape((-1, 2, 40, 2, 3))  # separate in to 2 regions
            val_x_grouped = val_x.reshape((-1, 2, 40, 2, 3))
            grouping_dict = {
                "all": {"swap": (1, 2), "shape": 12, "names": grouped_feature_names},
                # "hemi": {"swap": (3, 4), "shape": 2, "names": ["L", "R"]},
                # "scalar": {"swap": None, "shape": 3, "names": ["FA", "MD", "MK"]},
                # "region": {"swap": (1, 4), "shape": 2, "names": ["ANT", "POST"]},
            }
            for name, info in grouping_dict.items():
                if info["swap"] is not None:
                    this_shaps = np.swapaxes(shaps_grouped, info["swap"][0], info["swap"][1])
                    this_val = np.swapaxes(val_x_grouped, info["swap"][0], info["swap"][1])
                else:
                    this_shaps = shaps_grouped
                    this_val = val_x_grouped
                this_shaps = this_shaps.reshape((-1, info["shape"]))
                this_val = this_val.reshape((-1, info["shape"]))
                this_df = pd.DataFrame()
                for ii, feat_name in enumerate(grouped_feature_names):
                    section = feat_name.split("_")[2]
                    hemi = feat_name.split("_")[1]
                    scalar = feat_name.split("_")[0]
                    this_df = pd.concat((this_df, pd.DataFrame({
                        'Feature Value': this_val[:, ii],
                        'SHAP Value': this_shaps[:, ii],
                        "TP": scalar,
                        "Hemi": hemi,
                        "Section" : section})), ignore_index=True)
                bin_to_use = alt.Bin(extent=[-3.3, 3.3], step=0.3)
                shap_extent_to_use = alt.Scale(domain=(-0.03, 0.03))
                alt_shap_chart = alt.Chart().mark_line(size=line_size).encode(
                    x=alt.X("Feature Value", bin=bin_to_use),
                    y=alt.Y('mean(SHAP Value)', scale=shap_extent_to_use), color="TP")
                # alt_shap_chart_err = alt.Chart().mark_line(strokeDash=[1,1]).encode(
                #     x=alt.X("Feature Value", bin=bin_to_use),
                #     y='stdev(SHAP Value)', color="TP")
                alt_shap_chart_err = alt.Chart().mark_errorband(extent='iqr').encode(
                    x=alt.X("Feature Value", bin=bin_to_use),
                    y=alt.Y('SHAP Value', scale=shap_extent_to_use), color="TP"
                )
                alt_shap_chart = alt.layer(alt_shap_chart, alt_shap_chart_err, data=this_df).facet(
                    column=alt.Column('Hemi', header=alt.Header(titleFontSize=font_size, labelFontSize=font_size)),
                    row=alt.Row("Section", header=alt.Header(titleFontSize=font_size, labelFontSize=font_size)))
                alt_shap_chart.configure_axis(
                    labelFontSize=font_size,
                    titleFontSize=font_size
                    ).configure_legend(
                        labelFontSize=font_size,
                        titleFontSize=font_size,
                        labelLimit=0,
                        symbolStrokeWidth=10,
                    ).configure_title(
                        fontSize=font_size
                        ).save(
                            f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}_{name}_grouped.html")
                this_df.to_csv(f"output/pheno_match/{pheno_of_interest}/{use_model}_models/{b_name}_{name}_grouped.csv")

pheno_of_interest = "glauc_sec_nocut"  # can be: glauc, glauc_nom, glauc_sec, glauc_mega, logmar, ambly, alzh, any, age, glauc_sec_nocut, age_nocut
os.makedirs(f"output/pheno_match/{pheno_of_interest}/", exist_ok=True)
matcher = "mdm" # mdm, psm
global_recompute = False
use_model = "resnet"

pheno_dict = {
    "21003-2.0": "age",
    "31-0.0": "sex",
    "21000-0.0": "ethnicity",
    "3062-2.0": "FVC",
    "3064-2.0": "PEF",
    "189-0.0": "TDI"}
if pheno_of_interest == "glauc_nom":
    pheno_names = ["age", "sex", "ethnicity"]
elif pheno_of_interest == "glauc_sec":
    # pheno_names = list(pheno_dict.values())
    pheno_names = ["age", "sex", "ethnicity", "TDI"]
elif pheno_of_interest == "glauc_sec_nocut":
    pheno_names = ["age", "sex", "ethnicity", "TDI"]
elif pheno_of_interest == "glauc_mega":
    pheno_names = list(pheno_dict.values())
else:
    pheno_names = ["age", "sex", "ethnicity", "FVC", "PEF"]
all_pheno_names = list(pheno_dict.values())

if __name__ == "__main__":
    if True or global_recompute or not op.exists(f"output/pheno_match/{pheno_of_interest}/match.csv"):
        match_pipe()
    if True or global_recompute or not op.exists(f"output/pheno_match/{pheno_of_interest}/wide.csv"):
        wide_pipe()
    dis_df_wide = pd.read_csv(f"output/pheno_match/{pheno_of_interest}/wide.csv")
    if True or global_recompute or not op.exists(f"output/pheno_match/{pheno_of_interest}/ase_info.html"):
        ase_info_pipe()
    if True or global_recompute or not op.exists(f"output/pheno_match/{pheno_of_interest}/profiles.csv"):
        profiles_pipe()
    if True or global_recompute or not op.exists(f"output/aci_for_paper/profiles_w_aci_{pheno_of_interest}_match.csv"):
        aci_pipe()
    if True or global_recompute or not op.exists(f"output/pheno_match/{pheno_of_interest}/fracridge_roc.html"):
        fracridge_pipe()
    if True or global_recompute or not op.exists(f"output/pheno_match/{pheno_of_interest}/cnn.html"):
        cnn_pipe()
