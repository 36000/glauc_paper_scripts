import pandas as pd
import altair as alt

b_name = "OR"
name = "all"
font_size = 20
line_size = 5

glauc = pd.read_csv(f"output/pheno_match/glauc_sec/resnet_models/{b_name}_{name}_grouped.csv")
age = pd.read_csv(f"output/pheno_match/age/resnet_models/{b_name}_{name}_grouped.csv")

glauc["Phenotype"] = "Glaucoma"
age["Phenotype"] = "Age"

all_shaps = pd.concat((glauc, age))
all_shaps["Section"] = all_shaps["Section"] + ", " + all_shaps["Hemi"]
all_shaps["Hemisphere"] = all_shaps["Hemi"]
all_shaps["0_LINE"] = 0

bin_to_use = alt.Bin(extent=[-2.5, 2.5], step=0.5)
shap_extent_to_use = alt.Scale(domain=(-0.01, 0.01))
alt_shap_chart = alt.Chart().mark_line(size=line_size).encode(
    x=alt.X("Feature Value", bin=bin_to_use),
    y=alt.Y('mean(SHAP Value)', scale=shap_extent_to_use), color="Phenotype")
alt_shap_chart_err = alt.Chart().mark_errorband(extent='ci').encode( # extent iqr
    x=alt.X("Feature Value", bin=bin_to_use),
    y=alt.Y('SHAP Value', scale=shap_extent_to_use), color="Phenotype"
)
b_line_chart_y = alt.Chart().mark_rule(color='black').encode(y=alt.Y('0_LINE', scale=shap_extent_to_use))
b_line_chart_x = alt.Chart().mark_rule(color='black').encode(x='0_LINE')

print(all_shaps)
alt_shap_chart = alt.layer(b_line_chart_y + b_line_chart_x, alt_shap_chart + alt_shap_chart_err, data=all_shaps).facet(
    column=alt.Column('TP', header=alt.Header(titleFontSize=font_size, labelFontSize=font_size)))
   # row=alt.Row("Hemisphere", header=alt.Header(titleFontSize=font_size, labelFontSize=font_size)))
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
            f"output/second_paper_plots/shap_comp.html")
