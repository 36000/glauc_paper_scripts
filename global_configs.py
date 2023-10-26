font_size = 40
line_size = 5

bundle_names_formal = {
    "L_OR": "Left OR",
    "R_OR": "Right OR",
    "fov_L": "Left Fov. OR",
    "fov_R": "Right Fov. OR",
    "mac_L": "Left Mac. OR",
    "mac_R": "Right Mac. OR",
    "periph_L": "Left Per. OR",
    "periph_R": "Right Per. OR"
}

directions_formal = {
    "fov": "A→P",
    "periph": "A→P",
    "mac": "A→P",
    "CST": "I→S",
    "UNC": "P→A",
    "crossfp": "A→P",
    "crossfm": "A→P",
    "LR": "A→P",
    "LR_fov": "A→P",
    "LR_mac": "A→P",
    "LR_per": "A→P",
    "CST_L": "I→S",
    "CST_R": "I→S",
    "LR_CST": "I→S",
    "UNC_L": "P→A",
    "UNC_R": "P→A",
    "L_OR": "A→P",
    "R_OR": "A→P",
    "LR_UNC": "P→A",
    "L": "A→P",
    "R": "A→P"
}

column_names_formal = {
    "fov": "Fovea",
    "periph": "Periphery",
    "mac": "Macula",
    "CST": "CST",
    "UNC": "UNC",
    "crossfp": "Fov(+) v. Per(-)",
    "crossfm": "Fov(+) v. Mac(-)",
    "LR": "L(+) v. R(-)",
    "LR_fov": "L(+) v. R(-), Fov.",
    "LR_mac": "L(+) v. R(-), Mac.",
    "LR_per": "L(+) v. R(-), Per.",
    "CST_L": "Left CST",
    "CST_R": "Right CST",
    "LR_CST": "L(+) v. R(-), CST",
    "UNC_L": "Left UNC",
    "UNC_R": "Right UNC",
    "LR_UNC": "L(+) v. R(-), UNC",
    "L_OR": "Left OR",
    "R_OR": "Right OR",
    "L": "Left",
    "R": "Right"
}

side_assoc = {
    "L": "fov",
    "R": "periph"
}

get_column = {
    "fov": {
        "L": "fov_L",
        "R": "fov_R"
    },
    "periph": {
        "L": "periph_L",
        "R": "periph_R"
    },
    "mac": {
        "L": "mac_L",
        "R": "mac_R"
    },
    "CST": {
        "L": "CST_L",
        "R": "CST_R"
    },
    "UNC": {
        "L": "UNC_L",
        "R": "UNC_R"
    },
    "crossfp": {
        "L": "fov_L_periph_L_aci",
        "R": "fov_R_periph_R_aci"
    },
    "crossfm": {
        "L": "fov_L_mac_L_aci",
        "R": "fov_R_mac_R_aci"
    },
    "LR": {
        "L": "fov_L_fov_R_aci",
        "R": "periph_L_periph_R_aci"
    },
    "LR_fov": {
        "LvR": "fov_L_fov_R_aci",
    },
    "LR_per": {
        "LvR": "periph_L_periph_R_aci",
    },
    "LR_mac": {
        "LvR": "mac_L_mac_R_aci",
    },
    "CST_L": {"LvR": "CST_L", "L": "CST_L"},
    "CST_R": {"LvR": "CST_R", "R": "CST_R"},
    "LR_CST": {"LvR": "CST_L_CST_R_aci"},
    "UNC_L": {"LvR": "UNC_L", "L": "UNC_L"},
    "UNC_R": {"LvR": "UNC_R", "R": "UNC_R"},
    "LR_UNC": {"LvR": "UNC_L_UNC_R_aci"},
    "L_OR": {"LvR": "L_OR", "L": "L_OR"},
    "R_OR": {"LvR": "R_OR", "R": "R_OR"},
}

scalar_domains = {
    "fa": (0.4, 0.7),
    "md": (0.8, 1.2),
    "mk": (0.85, 1.15)
}
ukbb_scalar_domains = {
    "fa": (0.4, 0.6),
    "md": (0.9, 1.2),
    "mk": (0.7, 1.0)
}

cst_scalar_domains = {
    "fa": (0.4, 0.8),
    "md": (0.8, 1.0),
    "mk": (0.8, 1.3)
}

unc_scalar_domains = {
    "fa": (0.3, 0.5),
    "md": (0.8, 1.0),
    "mk": (0.5, 0.8)
}

aci_domain = (-0.25, 0.25)
aci_domain2 = (-0.1, 0.1)
match_aci_domain = (-1.0, 1.0)

column_count = {
    "ukbb": {"L": 7, "R": 7, "LvR": 4},
    "match": {"L": 2, "R": 2, "LvR": 4},
    "profs": {"L": 2, "R": 2, "LvR": 2},
    "hcp_roi": {"L": 2},
    "hcp_reco": {"L": 2},
    "hcp_boot": {"L": 2},
    "hcp_bootroi": {"L": 2}
}

binning_to_title = {
    "sex": "Sex",
    "acuity_l": "Left Acuity",
    "site": "Imaging Site",
    "glauc": "Glaucoma Status",
    "is_dis": "Glaucoma Status"
}

binning_to_order = {
    "sex": ["Female", "Male"],
    "acuity_l": ["-0.34,-0.15", "-0.15,-0.1", "-0.1,0.0", "0.0,0.1", "0.1,1.2"],
    "site": [11025.0, 11026.0, 11027.0],
    "glauc": [0, 1],
    "is_dis": [0, 1]
}

binning_to_color = {
    "sex": False,
    "acuity_l": False,
    "site": True,
    "glauc": True,
    "is_dis": True
}

order = ["Left", "Right", "Fovea", "Macula", "Periphery"]
age_bin_order = ["45-51", "52-55", "56-59", "60-63", "64-67", "68-71", "72-81"]
