"""Reference checks for the property calculations. Run with: pytest"""
import pytest

import til_chem as tc

# Expected values from RDKit, cross-checked against published values for these drugs
REFERENCE = {
    "aspirin": ("CC(=O)Oc1ccccc1C(=O)O", dict(formula="C9H8O4", mw=180.16, tpsa=63.60, hbd=1, hba=4, rotb=2)),
    "ibuprofen": ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", dict(formula="C13H18O2", mw=206.28, tpsa=37.30, hbd=1, hba=2, rotb=4)),
    "caffeine": ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", dict(formula="C8H10N4O2", mw=194.19, tpsa=61.82, hbd=0, hba=6, rotb=0)),
    "paracetamol": ("CC(=O)Nc1ccc(O)cc1", dict(formula="C8H9NO2", mw=151.16, tpsa=49.33, hbd=2, hba=3, rotb=1)),
}


@pytest.mark.parametrize("name", REFERENCE)
def test_reference_values(name):
    smiles, exp = REFERENCE[name]
    r = tc.analyse(smiles, name)
    assert r.ok
    assert r.formula == exp["formula"]
    assert r.mw == pytest.approx(exp["mw"], abs=0.01)
    assert r.tpsa == pytest.approx(exp["tpsa"], abs=0.01)
    assert (r.hbd, r.hba, r.rotb) == (exp["hbd"], exp["hba"], exp["rotb"])
    assert r.lipinski_pass and r.veber_pass


def test_canonical_smiles_is_stable():
    a = tc.analyse("OC(=O)c1ccccc1OC(C)=O")
    b = tc.analyse("CC(=O)Oc1ccccc1C(=O)O")
    assert a.canonical_smiles == b.canonical_smiles


def test_beyond_rule_of_five():
    r = tc.analyse("CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O", "atorvastatin")
    assert r.lipinski_violations == ["Molecular weight", "LogP"]
    assert not r.lipinski_pass


def test_invalid_and_empty_input():
    assert not tc.analyse("C1CC").ok          # unclosed ring
    assert not tc.analyse("C(C)(C)(C)(C)C").ok  # five-valent carbon
    assert not tc.analyse("   ").ok


def test_salt_flagged():
    assert tc.analyse("CC(=O)Oc1ccccc1C(=O)[O-].[Na+]").multiple_fragments


def test_batch_parsing():
    pairs = tc.parse_batch("# comment\nCCO ethanol\nWater,O\nc1ccccc1\n")
    assert pairs == [("CCO", "ethanol"), ("O", "Water"), ("c1ccccc1", "")]


def test_drawing_returns_png():
    png = tc.draw_png("c1ccccc1")
    assert png and png[:4] == b"\x89PNG"
