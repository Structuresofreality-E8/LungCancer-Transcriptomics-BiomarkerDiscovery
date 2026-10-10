from pathlib import Path
import streamlit as st
import pandas as pd
import json
ROOT = Path(__file__).resolve().parent
st.set_page_config(page_title="LUAD research explorer",layout="wide")
st.title("LUAD tissue biomarker research explorer")
st.caption("Retrospective research results; scores are not clinical probabilities.")
st.json(json.loads((ROOT/"results/metrics.json").read_text()))
sig = pd.read_csv(ROOT/"results/locked_signature.csv",index_col=0)
st.subheader("Locked signature")
st.dataframe(sig)
rep = pd.read_csv(ROOT/"results/external_paired_gene_replication.csv",index_col=0)
gene = st.selectbox("Inspect gene",sig.index.tolist())
st.dataframe(rep.loc[[gene]])
for name in ["external_confusion_matrix", "validation_curves", "locked_coefficients", "geo_signature_heatmap", "reactome_dotplot"]:
    p = ROOT/"figures"/(name+".png")
    if p.exists(): st.image(str(p))
st.markdown((ROOT/"LIMITATIONS.md").read_text())
