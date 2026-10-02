"""TIL Molecular Learning Workspace (v0.1)

A free, browser-based molecular workspace for members of The Insilico Lab.
Run locally with:  streamlit run app.py
"""
from __future__ import annotations

import io
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests
import streamlit as st

import til_chem as tc

APP_NAME = "TIL Molecular Learning Workspace"
VERSION = "0.1"
SITE_URL = "https://www.theinsilicolab.org/home"
DATA_DIR = Path(__file__).parent / "data"

st.set_page_config(page_title=APP_NAME, page_icon="🧪", layout="wide")


# ---------------------------------------------------------------- helpers
@st.cache_data(show_spinner=False)
def reference_molecules() -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / "reference_molecules.csv")


@st.cache_data(show_spinner=False, ttl=24 * 3600)
def pubchem_smiles(name: str) -> tuple[str | None, str]:
    """Look up a compound name on PubChem. Returns (smiles, message)."""
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/{}/JSON"
    for prop in ("SMILES", "IsomericSMILES"):
        try:
            r = requests.get(base.format(quote(name), prop), timeout=10)
        except requests.RequestException:
            return None, "PubChem could not be reached. Check your connection or paste the SMILES instead."
        if r.status_code == 404:
            return None, f"PubChem has no compound called “{name}”. Check the spelling or try another name."
        if r.ok:
            props = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
            smi = props.get("SMILES") or props.get("IsomericSMILES")
            if smi:
                return smi, f"Found on PubChem (CID {props.get('CID', '?')})."
    return None, "PubChem returned an unexpected answer. Paste the SMILES instead."


def results_csv(rows: list[dict]) -> bytes:
    return pd.DataFrame(rows).to_csv(index=False).encode("utf-8")


def rule_table(res: tc.MoleculeResult, limits: dict) -> pd.DataFrame:
    rows = []
    for label, (attr, limit) in limits.items():
        value = getattr(res, attr)
        unit = " Å²" if attr == "tpsa" else (" g/mol" if attr == "mw" else "")
        shown = f"{value:.2f}" if isinstance(value, float) else str(value)
        rows.append({
            "Rule": f"{label} ≤ {limit:g}{unit}",
            "This molecule": shown,
            "Result": "Pass" if value <= limit else "Fail",
        })
    return pd.DataFrame(rows)


if "saved" not in st.session_state:
    st.session_state.saved = []  # rows the trainee chose to keep this session
if "smiles" not in st.session_state:
    st.session_state.smiles = "CC(=O)Oc1ccccc1C(=O)O"
    st.session_state.mol_name = "Aspirin"


# ---------------------------------------------------------------- pages
def home():
    st.title(APP_NAME)
    st.caption(f"Version {VERSION}. Free for members of The Insilico Lab.")
    st.write(
        "Draw on the same tools used in computational drug discovery without installing anything. "
        "Type a molecule as SMILES or look it up by name, see its structure, calculate its properties, "
        "check it against Lipinski's rule of five, and download your results for assignments."
    )
    c1, c2, c3 = st.columns(3)
    with c1:
        st.subheader("Analyse one molecule")
        st.write("Structure, properties and drug-likeness rules for a single compound.")
        st.page_link(pages["workspace"], label="Open the workspace", icon="🧪")
    with c2:
        st.subheader("Compare many")
        st.write("Paste a list or upload a CSV and get one table for all of them.")
        st.page_link(pages["compare"], label="Compare molecules", icon="📊")
    with c3:
        st.subheader("Learn the terms")
        st.write("Short explanations of each property and rule, written for beginners.")
        st.page_link(pages["learn"], label="Read the guide", icon="📘")
    st.divider()
    st.write(f"Part of [The Insilico Lab]({SITE_URL}). Calculations use the open-source RDKit toolkit.")


def workspace():
    st.title("Molecule workspace")

    refs = reference_molecules()
    with st.container(border=True):
        tab_smiles, tab_name, tab_example = st.tabs(["Paste SMILES", "Search by name", "Use an example"])
        with tab_smiles:
            smi = st.text_input("SMILES", value=st.session_state.smiles, key="smiles_box",
                                help="A text code for a molecule, e.g. CC(=O)Oc1ccccc1C(=O)O for aspirin.")
            nm = st.text_input("Name (optional)", value=st.session_state.mol_name, key="name_box")
            if st.button("Analyse", type="primary"):
                st.session_state.smiles, st.session_state.mol_name = smi.strip(), nm.strip()
        with tab_name:
            query = st.text_input("Compound name", placeholder="e.g. imatinib, remdesivir, quercetin")
            if st.button("Search PubChem") and query.strip():
                with st.spinner("Searching PubChem…"):
                    found, msg = pubchem_smiles(query.strip())
                if found:
                    st.session_state.smiles, st.session_state.mol_name = found, query.strip().title()
                    st.success(msg)
                else:
                    st.error(msg)
        with tab_example:
            choice = st.selectbox("Example molecules", refs["name"], index=None, placeholder="Choose one")
            if choice:
                row = refs[refs["name"] == choice].iloc[0]
                if st.session_state.mol_name != choice:
                    st.session_state.smiles, st.session_state.mol_name = row["smiles"], choice
                    st.rerun()

    res = tc.analyse(st.session_state.smiles, st.session_state.mol_name)
    if not res.ok:
        st.error(res.error)
        return

    left, right = st.columns([2, 3], gap="large")
    with left:
        st.subheader(res.name or "Your molecule")
        png = tc.draw_png(res.canonical_smiles)
        if png:
            st.image(png, width="stretch")
        st.code(res.canonical_smiles, language=None)
        st.caption(f"Formula {res.formula}. Canonical SMILES as written by RDKit.")
        if res.multiple_fragments:
            st.warning("This SMILES has more than one fragment, often a salt or solvent. "
                       "Properties are calculated for everything together; remove the counter-ion to "
                       "describe the drug itself.")

    with right:
        st.subheader("Properties")
        m = st.columns(3)
        m[0].metric("Mol. weight (g/mol)", f"{res.mw:.2f}")
        m[1].metric("LogP (Crippen)", f"{res.logp:.2f}")
        m[2].metric("TPSA (Å²)", f"{res.tpsa:.1f}")
        m = st.columns(3)
        m[0].metric("H-bond donors", res.hbd)
        m[1].metric("H-bond acceptors", res.hba)
        m[2].metric("Rotatable bonds", res.rotb)
        m = st.columns(3)
        m[0].metric("Heavy atoms", res.heavy_atoms)
        m[1].metric("Aromatic rings", res.aromatic_rings)
        m[2].metric("Fraction sp3 C", f"{res.fsp3:.2f}")

        st.subheader("Lipinski's rule of five")
        n = len(res.lipinski_violations)
        if n == 0:
            st.success("Passes all four rules.")
        elif n == 1:
            st.info(f"One violation ({res.lipinski_violations[0]}). One is usually allowed, so this still counts as drug-like.")
        else:
            st.error(f"{n} violations: {', '.join(res.lipinski_violations)}. Oral absorption may be poor, "
                     "though some approved drugs, such as many antibiotics and cyclosporin, sit outside these rules.")
        st.dataframe(rule_table(res, tc.LIPINSKI_LIMITS), hide_index=True, width="stretch")

        with st.expander("Veber's rules for oral bioavailability"):
            st.dataframe(rule_table(res, tc.VEBER_LIMITS), hide_index=True, width="stretch")

    st.divider()
    row = res.as_row()
    c1, c2, c3 = st.columns([1, 1, 2])
    if c1.button("Add to my results"):
        if any(r["Canonical SMILES"] == row["Canonical SMILES"] for r in st.session_state.saved):
            st.toast("Already in your results.")
        else:
            st.session_state.saved.append(row)
            st.toast(f"Added {res.name or 'molecule'} to your results.")
    c2.download_button("Download this molecule (CSV)", results_csv([row]),
                       file_name=f"{(res.name or 'molecule').replace(' ', '_')}_properties.csv", mime="text/csv")

    if st.session_state.saved:
        st.subheader(f"My results ({len(st.session_state.saved)})")
        st.dataframe(pd.DataFrame(st.session_state.saved), hide_index=True, width="stretch")
        d1, d2, _ = st.columns([1, 1, 2])
        d1.download_button("Download all results (CSV)", results_csv(st.session_state.saved),
                           file_name="til_results.csv", mime="text/csv", type="primary")
        if d2.button("Clear results"):
            st.session_state.saved = []
            st.rerun()
        st.caption("Results are kept only while this page is open. Download them before you leave.")


def compare():
    st.title("Compare molecules")
    st.write("Put one molecule per line as `SMILES name`, or `name,SMILES`. You can also upload a CSV "
             "with columns called `name` and `smiles`.")
    sample = "\n".join(f"{r.smiles} {r.name}" for r in reference_molecules().head(5).itertuples())
    text = st.text_area("Molecules", value=sample, height=180)
    upload = st.file_uploader("Or upload a CSV", type=["csv"])

    pairs: list[tuple[str, str]] = []
    if upload is not None:
        try:
            df = pd.read_csv(upload)
            cols = {c.lower().strip(): c for c in df.columns}
            if "smiles" not in cols:
                st.error("The CSV needs a column called “smiles”.")
                return
            names = df[cols["name"]].astype(str) if "name" in cols else [""] * len(df)
            pairs = list(zip(df[cols["smiles"]].astype(str), names))
        except Exception as exc:  # noqa: BLE001
            st.error(f"That file could not be read as a CSV ({exc}).")
            return
    else:
        pairs = tc.parse_batch(text)

    if not pairs:
        st.info("Add at least one molecule above.")
        return
    if len(pairs) > 500:
        st.warning("Only the first 500 molecules are shown. Use the Colab notebooks for larger sets.")
        pairs = pairs[:500]

    results = [tc.analyse(s, n) for s, n in pairs]
    rows = [r.as_row() for r in results]
    bad = [r for r in results if not r.ok]
    good = pd.DataFrame([r.as_row() for r in results if r.ok])

    if bad:
        st.warning(f"{len(bad)} of {len(results)} lines could not be read: "
                   + ", ".join(b.name or b.input_smiles[:25] for b in bad[:6]) + ("…" if len(bad) > 6 else ""))
    if good.empty:
        return

    passed = (good["Lipinski"] == "Pass").sum()
    m = st.columns(3)
    m[0].metric("Molecules", len(good))
    m[1].metric("Pass Lipinski", f"{passed} of {len(good)}")
    m[2].metric("Median MW", f"{good['MW (g/mol)'].median():.1f}")

    front = ["Name", "MW (g/mol)", "LogP (Crippen)", "TPSA (Å²)", "HBD", "HBA", "Rotatable bonds",
             "Lipinski violations", "Lipinski", "Veber"]
    shown = good[front + [c for c in good.columns if c not in front and c != "Error"]]
    st.dataframe(shown, hide_index=True, width="stretch")
    st.download_button("Download table (CSV)", results_csv(rows), file_name="til_comparison.csv",
                       mime="text/csv", type="primary")

    st.subheader("Molecular weight against LogP")
    st.caption("The rule-of-five region is MW ≤ 500 and LogP ≤ 5, the lower-left of the chart.")
    chart = good.assign(Lipinski=good["Lipinski"])
    st.scatter_chart(chart, x="LogP (Crippen)", y="MW (g/mol)", color="Lipinski", height=380)


def learn():
    st.title("What the numbers mean")
    st.write("A quick guide to every value in the workspace. These are rules of thumb for oral drugs, "
             "not laws: they help you prioritise molecules, not reject them outright.")
    guide = [
        ("SMILES", "A line of text that describes a molecule's atoms and bonds. Lowercase letters are aromatic "
         "atoms, numbers open and close rings, brackets show branches. RDKit rewrites your input into a "
         "canonical SMILES so the same molecule always gets the same text."),
        ("Molecular weight", "The mass of one mole of the molecule in g/mol. Larger molecules generally cross "
         "membranes less easily. Lipinski's limit is 500."),
        ("LogP", "How much the molecule prefers oil (octanol) over water, on a log scale. Higher means more "
         "lipophilic. This app uses the Crippen method; other tools such as PubChem's XLogP3 or SwissADME "
         "give somewhat different values, so always say which method you used. Lipinski's limit is 5."),
        ("TPSA", "Topological polar surface area: the surface belonging to polar atoms (mostly N and O and "
         "their hydrogens), in Å². High TPSA tends to mean poor membrane permeability. Veber's limit is 140 Å²; "
         "for crossing the blood–brain barrier, values under about 90 Å² are usually wanted."),
        ("H-bond donors and acceptors", "Counted here with Lipinski's original definitions: donors are NH and "
         "OH groups, acceptors are all N and O atoms. Other software counts them differently, which is why "
         "numbers can disagree between tools. Limits: 5 donors, 10 acceptors."),
        ("Rotatable bonds", "Single bonds outside rings that let the molecule change shape. Very flexible "
         "molecules tend to have poorer oral bioavailability. Veber's limit is 10. In docking, more rotatable "
         "bonds also make the search harder."),
        ("Fraction sp3 C", "The share of carbon atoms that are sp3 (saturated). Flat, aromatic-heavy molecules "
         "have low values; a moderate value is often linked with better solubility."),
        ("Lipinski's rule of five", "Most orally active drugs have no more than one violation of: MW ≤ 500, "
         "LogP ≤ 5, donors ≤ 5, acceptors ≤ 10. Natural products, antibiotics and some newer drug classes "
         "deliberately sit beyond these rules, so treat a failure as a flag to look closer."),
    ]
    for term, text in guide:
        with st.expander(term):
            st.write(text)
    st.caption("References: Lipinski et al., Adv. Drug Deliv. Rev. 1997; Veber et al., J. Med. Chem. 2002; "
               "Wildman and Crippen, J. Chem. Inf. Comput. Sci. 1999; Ertl et al., J. Med. Chem. 2000.")


# ---------------------------------------------------------------- navigation
pages = {
    "home": st.Page(home, title="Home", icon="🏠", default=True),
    "workspace": st.Page(workspace, title="Molecule workspace", icon="🧪", url_path="workspace"),
    "compare": st.Page(compare, title="Compare molecules", icon="📊", url_path="compare"),
    "learn": st.Page(learn, title="What the numbers mean", icon="📘", url_path="learn"),
}
nav = st.navigation(list(pages.values()))
with st.sidebar:
    st.markdown(f"**The Insilico Lab**  \n[theinsilicolab.org]({SITE_URL})")
    st.caption(f"{APP_NAME} v{VERSION}")
nav.run()
