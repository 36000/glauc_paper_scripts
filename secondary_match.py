from primary_match import auc_score, make_roc_plot, print_delong_table, get_auc_unc, formal_or_names
import tensorflow as tf
import pandas as pd
import numpy as np
from scipy.optimize import linear_sum_assignment
from tqdm import tqdm
import os.path as op
import os
from sklearn.metrics import roc_auc_score, roc_curve
import sys

pheno_dict = {
    "21003-2.0": "age",
    "31-0.0": "sex",
    "21000-0.0": "ethnicity",
    "3062-2.0": "FVC",
    "3064-2.0": "PEF",
    "189-0.0": "TDI"}
pheno_names = ["age", "sex", "ethnicity", "TDI"]


if len(sys.argv) > 1:
    pheno_of_interest = sys.argv[1]
else:
    pheno_of_interest = "amd_rob" # amd_rob, tenano, glauc_rob, amd_rob_age

if len(sys.argv) > 2:
    from_folder_suffix = sys.argv[2]
else:
    from_folder_suffix = "" # "", "_nocut"


if pheno_of_interest == "glauc_rob" or pheno_of_interest == "amd_rob_age":
    from_folder = "age"
else:
    from_folder = "glauc_sec"

from_folder = from_folder + from_folder_suffix
pheno_of_interest = pheno_of_interest + from_folder_suffix

pheno = pd.read_csv("output/pheno.csv", low_memory=False)
os.makedirs(f"output/pheno_match/{pheno_of_interest}/", exist_ok=True)

if "amd_rob" in pheno_of_interest:
    dis_subs = pheno[
        (pheno["6148-2.0"] == 5)
        | (pheno["6148-3.0"] == 5)].eid.unique()

    # exclude acuity worse than 0.3
    if from_folder_suffix != "_nocut":
        for acuity_f in ["5201-0.0", "5201-1.0", "5208-0.0", "5208-1.0"]:
            pheno = pheno[(pheno[acuity_f] <= 0.3) | np.isnan(pheno[acuity_f])]

    # exclude other retinal disease
    pheno = pheno[
        (pheno["6148-2.0"] == -7)
        | (pheno["6148-3.0"] == -7)
        | pheno.eid.isin(dis_subs)]

    # exclude used controls
    pheno = pheno[~pheno.eid.isin(np.loadtxt(f"output/pheno_match/{from_folder}/ctrl_sub.txt"))]

    pheno = pheno.rename(columns=pheno_dict)
    pheno = pheno[["eid", *pheno_names]]

    test = pheno[pheno.eid.isin(dis_subs)]
    control = pheno[~pheno.eid.isin(dis_subs)]

    test.loc[:, f'{pheno_of_interest}_status'] = 1
    control.loc[:, f'{pheno_of_interest}_status'] = 0
    threshold = 0.3

elif "tenano" in pheno_of_interest:
    # exclude acuity worse than 0.3
    for acuity_f in ["5201-0.0", "5201-1.0", "5208-0.0", "5208-1.0"]:
        pheno = pheno[(pheno[acuity_f] <= 0.3) | np.isnan(pheno[acuity_f])]

    # exclude other retinal disease
    pheno = pheno[
        (pheno["6148-2.0"] == -7)
        | (pheno["6148-3.0"] == -7)]

    # exclude used controls
    pheno = pheno[~pheno.eid.isin(np.loadtxt(f"output/pheno_match/{from_folder}/ctrl_sub.txt"))]

    pheno = pheno.rename(columns=pheno_dict)
    pheno = pheno[["eid", *pheno_names]]

    test = pheno.copy()
    test = test[test["age"] == 70] # reduce linear_sum_assignment time
    control = pheno.copy()
    control["age"] = control["age"] + 10

    test.loc[:, f'{pheno_of_interest}_status'] = 1
    control.loc[:, f'{pheno_of_interest}_status'] = 0
    threshold = 0.3

elif "glauc_rob" in pheno_of_interest:
    dis_subs = pheno[np.logical_not(
        np.isnan(pheno["6119-0.0"]) &
        np.isnan(pheno["6119-1.0"]) &
        np.isnan(pheno["6119-2.0"]) &
        np.isnan(pheno["6119-3.0"]))]
    dis_subs = dis_subs[np.logical_and(dis_subs["21003-2.0"] >= 60, dis_subs["21003-2.0"] < 65)]
    dis_subs = dis_subs.eid.unique()
    #print(len(pheno))
    # exclude acuity worse than 0.3
    for acuity_f in ["5201-0.0", "5201-1.0", "5208-0.0", "5208-1.0"]:
        pheno = pheno[(pheno[acuity_f] <= 0.3) | np.isnan(pheno[acuity_f])]

    # exclude other retinal disease
    pheno = pheno[
        (pheno["6148-2.0"] == -7)
        | (pheno["6148-3.0"] == -7)
        | pheno.eid.isin(dis_subs)]

    # exclude glauc
    pheno = pheno[
        (np.isnan(pheno["6119-0.0"]) &
         np.isnan(pheno["6119-1.0"]) &
         np.isnan(pheno["6119-2.0"]) &
         np.isnan(pheno["6119-3.0"])) | pheno.eid.isin(dis_subs)]

    # exclude used controls
    pheno = pheno[~pheno.eid.isin(np.loadtxt(f"output/pheno_match/{from_folder}/ctrl_sub.txt"))]

    pheno = pheno.rename(columns=pheno_dict)
    pheno = pheno[["eid", *pheno_names]]

    test = pheno[pheno.eid.isin(dis_subs)]
    control = pheno[~pheno.eid.isin(dis_subs)]

    test.loc[:, f'{pheno_of_interest}_status'] = 1
    control.loc[:, f'{pheno_of_interest}_status'] = 0
    threshold = 0.3
else:
    raise ValueError("bad pheno of interest")


print(test)
print(control)


test_wo_eid = test.drop(["eid", f'{pheno_of_interest}_status'], axis=1).to_numpy()
print("Original test len: " + str(test_wo_eid.shape[0]))
ctrl_wo_eid = control.drop(["eid", f'{pheno_of_interest}_status'], axis=1).to_numpy()
print("Original ctrl len: " + str(ctrl_wo_eid.shape[0]))
test_wo_eid[np.isnan(test_wo_eid)] = -1
ctrl_wo_eid[np.isnan(ctrl_wo_eid)] = -1

v_inv = np.linalg.inv(np.cov(np.concatenate((test_wo_eid, ctrl_wo_eid), axis=0).T, ddof=0))
diff = test_wo_eid[:, np.newaxis] - ctrl_wo_eid[np.newaxis, :]
nbrs = np.sqrt(np.sum((diff@v_inv)*diff, axis=2))

# nbrs2 = np.zeros((test_wo_eid.shape[0], ctrl_wo_eid.shape[0]))
# for i in range(test_wo_eid.shape[0]):
#     for j in range(ctrl_wo_eid.shape[0]):
#         nbrs2[i, j] = mahalanobis(test_wo_eid[i], ctrl_wo_eid[j], v_inv)

row_ind, col_ind = linear_sum_assignment(nbrs)

filtered_row_ind = [] # remove matches that are too bad
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
for eids, frame in zip([test_eid, ctrl_eid], [test, control]):
    for matchid, eid in enumerate(eids):
        _dict = frame[frame.eid==eid].iloc[0].to_dict()
        _dict["match_id"] = matchid
        match_df = match_df.append(_dict, ignore_index=True)
for intcolumn in ["match_id", "eid", f'{pheno_of_interest}_status']:
    match_df[intcolumn] = match_df[intcolumn].astype(int)
print("Len of matched eid: " + str(len(test_eid)))
match_df.to_csv(f"output/pheno_match/{pheno_of_interest}/match.csv", index=False)

dis_df = pd.read_csv(f"output/pheno_match/{pheno_of_interest}/match.csv")
dis_df = dis_df.drop(columns=["scores", "record_id"], errors="ignore")

dis_df_wide = dis_df.pivot(
    index="match_id", columns=f'{pheno_of_interest}_status',
    values=[*pheno_names, "eid"]).reset_index()
dis_df_wide.columns = dis_df_wide.columns.map(str).map(''.join).astype(str).str.replace("'", "")
print(dis_df_wide)
dis_df_wide = dis_df_wide.drop_duplicates(["(eid, 0)"])
print(dis_df_wide)
dis_df_wide.to_csv(f"output/pheno_match/{pheno_of_interest}/wide.csv", index=False)

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
print(len(matched_profiles.subjectID_ctrl.unique()))
print(len(matched_profiles))
matched_profiles.to_csv(f"output/pheno_match/{pheno_of_interest}/profiles.csv")

dataframe = pd.read_csv(f"output/pheno_match/{pheno_of_interest}/profiles.csv")
all_subjects = dataframe.subjectID_ctrl.unique()
bd = {
    "OR": ["L_OR",  "R_OR"],
    "CST": ["CST_L", "CST_R"],
    "UNC": ["UNC_L", "UNC_R"]}

Xs = {}
already_exists = True
for b_name, bundles in bd.items():
    Xs[b_name] = {}
    if not op.exists((
            f"output/pheno_match/{pheno_of_interest}/"
            f"X_{b_name}.npy")):
        already_exists = False

if True or not already_exists:
    for b_name, bundles in bd.items():
        Xs[b_name] = np.zeros((len(all_subjects), 2, len(bundles), 3, 100))
    for ii, subjectID in enumerate(tqdm(all_subjects, leave=False)):
        sub_df = dataframe[dataframe.subjectID_ctrl == subjectID]
        for nodeID in tqdm(range(100), leave=False):
            node_df = sub_df[sub_df.nodeID == nodeID]
            for b_name, bundles in tqdm(bd.items(), leave=False):
                for ll, bundle in enumerate(tqdm(bundles, leave=False)):
                    for jj, scalar in enumerate(tqdm(["dki_fa", "dki_md", "dki_mk"], leave=False)):
                        for kk, suf in enumerate(tqdm(["ctrl", "dis"], leave=False)):
                            Xs[b_name][ii, kk, ll, jj, nodeID] = node_df[
                                f"{scalar}_{bundle}_{suf}"].to_numpy()[0]
    for b_name in bd.keys():
        np.save((
            f"output/pheno_match/{pheno_of_interest}/"
            f"X_{b_name}.npy"), Xs[b_name])

datasets = {}
raw_dataset = {}
batch_size = 128

for b_name in bd.keys():
    datasets[b_name] = {}
    raw_dataset[b_name] = {}
    this_X = np.load(
        f"output/pheno_match/{pheno_of_interest}/X_{b_name}.npy")

    target = np.zeros((this_X.shape[0], 2))
    target[:, 1] = 1
    target = target.flatten()

    this_X[
        np.any(np.isnan(this_X), axis=(2, 3, 4))] =\
        np.nanmean(this_X, axis=(0, 1))

    this_X = this_X[:, :, :, :, 10:90]
    this_X = this_X.reshape((-1, 6, 80))
    this_X = np.swapaxes(this_X, 1, 2)

    raw_dataset[b_name] = this_X

    this_std = np.std(this_X, axis=0)
    this_mean = np.mean(this_X, axis=0)
    this_X = (this_X-this_mean)/this_std

    dataset = tf.data.Dataset.from_tensor_slices((this_X, target))
    dataset = dataset.batch(batch_size)
    datasets[b_name] = dataset

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
        f"output/pheno_match/{from_folder}/resnet_models/{b_name}.h5")
    target = np.concatenate([y for x, y in datasets[b_name]], axis=0)
    target_hat = np.asarray(model_cnn.predict(datasets[b_name])[:, 0])
    print(f"{b_name}: {auc_score(target, target_hat)}")
    fpr, tpr, tresh = roc_curve(target, target_hat)
    this_roc_df = pd.DataFrame()
    this_roc_df['FPR'] = fpr
    this_roc_df['TPR'] = tpr
    this_roc_df['dotted_line'] = np.asarray(tpr > 0.5).astype(float)
    auc, ci = get_auc_unc(target, target_hat)
    bundle_name_w_roc = formal_or_names[b_name] + f", AUC: {round(auc, 2)}"
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

make_roc_plot(
    f"{pheno_of_interest}/cnn_roc",
    roc_df, auc_df, b_order, b_order_auc, '1D Convolutional Resnet ROC')