# TIL Molecular Learning Workspace

Version 0.1 of the TIL CADD Learning Hub: a free, browser-based molecular workspace for members of [The Insilico Lab](https://www.theinsilicolab.org/home). Trainees can calculate molecular properties and check drug-likeness without installing Python or RDKit.

## What it does

- Accepts a molecule as SMILES, by name (PubChem lookup) or from a set of examples
- Draws the 2D structure
- Calculates molecular weight, LogP (Crippen), TPSA, H-bond donors and acceptors, rotatable bonds, heavy atoms, rings and fraction sp3 carbon
- Checks Lipinski's rule of five and Veber's rules
- Compares many molecules at once from a pasted list or CSV upload
- Exports results as CSV for assignments
- Explains every property in plain language

## Project layout

```
app.py                        Streamlit interface (pages and layout)
til_chem.py                   All RDKit calculations, kept separate so they can be tested and reused in Colab
data/reference_molecules.csv  Example molecules used in the app and in exercises
tests/test_til_chem.py        Reference-value tests (run with pytest)
.streamlit/config.toml        Theme and upload limits
requirements.txt              Python packages for deployment
```

## Run it on your own computer

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt pytest
pytest                           # all tests should pass
streamlit run app.py
```

## Put it online (Streamlit Community Cloud, free)

1. Create a GitHub repository (for example `theinsilicolab/til-workspace`) and upload these files.
2. Sign in at https://share.streamlit.io with the same GitHub account.
3. Choose **Create app**, select the repository, branch `main` and main file `app.py`.
4. Under **App URL**, pick a readable address such as `til-workspace.streamlit.app`.
5. Deploy. The first build takes a few minutes because RDKit is large.

## Connect it to theinsilicolab.org (Google Sites)

1. Open the site in the Google Sites editor and add a page, for example **Tools**.
2. Choose **Insert → Embed → By URL** and paste the app address with `?embed=true` at the end:
   `https://til-workspace.streamlit.app/?embed=true`
3. Drag the embed box to at least 1,200 px tall.
4. Also add a button that opens the app in its own tab. Phones and slow connections work better that way.

## Known limits of the free tier

- Apps on Streamlit Community Cloud go to sleep after a period with no visitors, and the first visitor waits while it wakes up. Open the app yourself a few minutes before a live session.
- The free tier has limited memory shared by everyone using the app at once. It is fine for small groups; test it before a full cohort works in it simultaneously.
- Results are kept only while a page is open. Nothing is stored on the server and no personal data is collected.

## Methods

Lipinski et al., Adv. Drug Deliv. Rev. 1997. Veber et al., J. Med. Chem. 2002. Wildman and Crippen, J. Chem. Inf. Comput. Sci. 1999 (LogP). Ertl et al., J. Med. Chem. 2000 (TPSA). RDKit, https://www.rdkit.org.
