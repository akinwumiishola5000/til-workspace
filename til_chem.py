"""Molecular property calculations for the TIL Molecular Learning Workspace.

Every calculation lives here, separate from the interface, so that trainees
can read it, reuse it in Colab, and so it can be tested on its own.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")  # keep RDKit parse warnings out of the app logs

# Lipinski's rule of five (Lipinski et al., Adv. Drug Deliv. Rev. 1997)
LIPINSKI_LIMITS = {
    "Molecular weight": ("mw", 500.0),
    "LogP": ("logp", 5.0),
    "H-bond donors": ("hbd", 5),
    "H-bond acceptors": ("hba", 10),
}
# Veber's oral bioavailability rules (Veber et al., J. Med. Chem. 2002)
VEBER_LIMITS = {
    "Rotatable bonds": ("rotb", 10),
    "TPSA": ("tpsa", 140.0),
}


@dataclass
class MoleculeResult:
    input_smiles: str
    name: str = ""
    ok: bool = False
    error: str = ""
    canonical_smiles: str = ""
    formula: str = ""
    mw: float = 0.0
    logp: float = 0.0
    tpsa: float = 0.0
    hbd: int = 0
    hba: int = 0
    rotb: int = 0
    heavy_atoms: int = 0
    rings: int = 0
    aromatic_rings: int = 0
    fsp3: float = 0.0
    lipinski_violations: list[str] = field(default_factory=list)
    veber_violations: list[str] = field(default_factory=list)
    multiple_fragments: bool = False

    @property
    def lipinski_pass(self) -> bool:
        # The usual convention: no more than one violation is allowed
        return len(self.lipinski_violations) <= 1

    @property
    def veber_pass(self) -> bool:
        return not self.veber_violations

    def as_row(self) -> dict:
        if not self.ok:
            return {"Name": self.name, "Input SMILES": self.input_smiles, "Error": self.error}
        return {
            "Name": self.name,
            "Input SMILES": self.input_smiles,
            "Canonical SMILES": self.canonical_smiles,
            "Formula": self.formula,
            "MW (g/mol)": round(self.mw, 2),
            "LogP (Crippen)": round(self.logp, 2),
            "TPSA (Å²)": round(self.tpsa, 2),
            "HBD": self.hbd,
            "HBA": self.hba,
            "Rotatable bonds": self.rotb,
            "Heavy atoms": self.heavy_atoms,
            "Rings": self.rings,
            "Aromatic rings": self.aromatic_rings,
            "Fraction sp3 C": round(self.fsp3, 2),
            "Lipinski violations": len(self.lipinski_violations),
            "Lipinski": "Pass" if self.lipinski_pass else "Fail",
            "Veber": "Pass" if self.veber_pass else "Fail",
            "Error": "",
        }


def parse(smiles: str) -> Chem.Mol | None:
    smiles = (smiles or "").strip()
    if not smiles:
        return None
    return Chem.MolFromSmiles(smiles)


def analyse(smiles: str, name: str = "") -> MoleculeResult:
    res = MoleculeResult(input_smiles=(smiles or "").strip(), name=name.strip())
    if not res.input_smiles:
        res.error = "No SMILES given."
        return res
    mol = parse(res.input_smiles)
    if mol is None:
        res.error = "RDKit could not read this SMILES. Check brackets, ring numbers and atom valences."
        return res

    res.ok = True
    res.canonical_smiles = Chem.MolToSmiles(mol)
    res.formula = rdMolDescriptors.CalcMolFormula(mol)
    res.mw = Descriptors.MolWt(mol)
    res.logp = Crippen.MolLogP(mol)
    res.tpsa = rdMolDescriptors.CalcTPSA(mol)
    # Lipinski's original definitions: donors = NH + OH, acceptors = N + O
    res.hbd = rdMolDescriptors.CalcNumLipinskiHBD(mol)
    res.hba = rdMolDescriptors.CalcNumLipinskiHBA(mol)
    res.rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    res.heavy_atoms = mol.GetNumHeavyAtoms()
    res.rings = rdMolDescriptors.CalcNumRings(mol)
    res.aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    res.fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    res.multiple_fragments = len(Chem.GetMolFrags(mol)) > 1

    for label, (attr, limit) in LIPINSKI_LIMITS.items():
        if getattr(res, attr) > limit:
            res.lipinski_violations.append(label)
    for label, (attr, limit) in VEBER_LIMITS.items():
        if getattr(res, attr) > limit:
            res.veber_violations.append(label)
    return res


def draw_png(smiles: str, width: int = 420, height: int = 320) -> bytes | None:
    mol = parse(smiles)
    if mol is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.padding = 0.08
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def parse_batch(text: str) -> list[tuple[str, str]]:
    """Read one molecule per line: 'SMILES' or 'SMILES name' or 'name,SMILES'."""
    out = []
    for raw in (text or "").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "," in line:
            a, b = [p.strip() for p in line.split(",", 1)]
            # whichever side parses as a molecule is the SMILES
            if parse(a) is not None:
                out.append((a, b))
            else:
                out.append((b, a))
        else:
            parts = line.split(None, 1)
            out.append((parts[0], parts[1] if len(parts) > 1 else ""))
    return out
