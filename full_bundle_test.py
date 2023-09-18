import pandas as pd
import numpy as np
from numpy.random import default_rng
import os.path as op
import tensorflow as tf
from primary_match import auc_score, formal_or_names, print_delong_table, make_roc_plot
from sklearn.metrics import roc_auc_score, roc_curve

pheno_of_interest = "glauc_sec"
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
    # "train": all_subjects[training_idx],
    "test": all_subjects[test_idx],}
    # "val": all_subjects[val_idx]}

bd = {
    # "FOV": ["fov_L",  "fov_R"],
    # "MAC": ["mac_L",  "mac_R"],
    # "PERIPH": ["periph_L", "periph_R"],
    "OR": ["L_OR", "R_OR"],
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

for b_name in bd.keys():
    for sett in subs.keys():
        Xs[b_name][sett] = np.load((
            f"output/pheno_match/{pheno_of_interest}/"
            f"fracridge_derivs/X_{b_name}_{sett}.npy"))

# remove nans
# c_valid = 0
idxs = None
for sett in subs.keys():
    for b_name, bundles in bd.items():
        other_new_idxs = np.all(np.equal(
           np.isnan(Xs[b_name][sett][:, 0, ...]),
           np.isnan(Xs[b_name][sett][:, 1, ...])), axis=(1, 2, 3))
        new_idxs = ~np.any(np.isnan(Xs[b_name][sett]), axis=(1, 2, 3, 4))
        if idxs is None:
            idxs = other_new_idxs
        else:
            idxs = np.logical_and(idxs, other_new_idxs)
    for b_name, bundles in bd.items():
        print(Xs[b_name][sett].shape)
        Xs[b_name][sett] = Xs[b_name][sett][idxs]
        Xs[b_name][sett][
            np.any(np.isnan(Xs[b_name][sett]), axis=(2, 3, 4))] =\
            np.nanmean(Xs[b_name][sett], axis=(0, 1))
        print(Xs[b_name][sett].shape)

        # Xs[b_name][sett][
        #     np.any(np.isnan(Xs[b_name][sett]), axis=(2, 3, 4))] =\
        #     np.nanmean(Xs[b_name][sett], axis=(0, 1))

datasets = {}
batch_size = 128

for b_name in bd.keys():
    datasets[b_name] = {}
    for sett in subs:
        this_X = Xs[b_name][sett]

        target = np.zeros((this_X.shape[0], 2))
        target[:, 1] = 1
        target = target.flatten()

        this_X = this_X[:, :, :, :, 10:90]
        this_X = this_X.reshape((-1, 6, 80))
        this_X = np.swapaxes(this_X, 1, 2)

        this_std = np.std(this_X, axis=0)
        this_mean = np.mean(this_X, axis=0)
        this_X = (this_X-this_mean)/this_std

        dataset = tf.data.Dataset.from_tensor_slices((this_X, target))
        dataset = dataset.batch(batch_size)
        datasets[b_name][sett] = dataset

roc_df = pd.DataFrame()
auc_df = pd.DataFrame()
b_order = []
b_order_auc = []
hats = []
for b_name in bd.keys():
    model_cnn = tf.keras.models.load_model(
        f"output/pheno_match/{pheno_of_interest}/resnet_models/{b_name}.h5")
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

    # tf.keras.utils.plot_model(
    #     model_cnn,
    #     to_file=f"output/pheno_match/{pheno_of_interest}/resnet_models/{b_name}_model.png",
    #     show_layer_names=False,
    #     show_shapes=False)

print(target.shape)
print(hats[0].shape)
print_delong_table(target, hats)

make_roc_plot(
    f"{pheno_of_interest}/cnn_roc_only_full",
    roc_df, auc_df, b_order, b_order_auc, '1D Convolutional Resnet ROC')