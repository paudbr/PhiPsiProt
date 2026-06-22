#!/usr/bin/env python3
"""
structural_metrics.py — calcula 6 métricas estructurales sobre un PDB de modelo
(de ESMFold o Chai-1) y escribe un TSV largo: id <TAB> metric <TAB> value.

Métricas:
  rosetta_dg    ΔG absoluto (REU) del modelo, ref2015, tras FastRelax restringido
                con restricciones de coordenadas sobre CA (coord_cst). lower=mejor.
  tm_score      TM-score del modelo vs PDB de referencia (input_pdb). higher=mejor.
  mean_plddt    media del B-factor del PDB (pLDDT por residuo). higher=mejor.
  n_cysteines   nº de residuos CYS. (normalmente gate, no peso)
  net_charge    carga neta a pH dado, pKa por entorno con PROPKA.
  net_charge_abs |net_charge| (para ponderar "cerca de neutro"). lower=mejor.
  sap_score     Spatial Aggregation Propensity (proxy de solubilidad). lower=mejor.

Uso:
  structural_metrics.py --id <cand_id> --model model.pdb [--ref ref.pdb]
                        [--ph 7.0] [--relax restrained|none|full]
                        [--out <id>_metrics.tsv]

Sin --ref, tm_score se omite (NA).
"""
import argparse
import sys
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np


# ───────────────────────── parser PDB/CIF ─────────────────────────────────────
def get_parser(path):
    """Devuelve el parser de Biopython adecuado según la extensión."""
    from Bio.PDB import PDBParser, MMCIFParser
    if str(path).lower().endswith(".cif"):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


# Aminoácidos estándar (3 letras). Para filtrar HETATM/residuos no proteicos.
_STD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    # variantes de protonación de histidina que algunos PDB usan
    "HID","HIE","HIP","MSE",
}


def _strip_hetatm(path):
    """Crea un PDB temporal SOLO con átomos de proteína (líneas ATOM de residuos
    estándar). Elimina HETATM, aguas, iones (ZN), ligandos (AZM), glicerol (GOL),
    etc. Devuelve la ruta del PDB limpio. Si la entrada es CIF, la parsea con
    Biopython y reescribe solo la proteína.

    Imprescindible para que el ΔG (y por tanto el ΔΔG) sea comparable entre el
    original cristalográfico (con cofactores) y los modelos (proteína sola).
    """
    import os
    out = "_clean_" + os.path.basename(str(path)).rsplit(".", 1)[0] + ".pdb"

    if str(path).lower().endswith(".cif"):
        # CIF: parsear y reescribir solo residuos estándar
        from Bio.PDB import PDBIO, Select
        parser = get_parser(path)
        st = parser.get_structure("s", path)

        class ProtSelect(Select):
            def accept_residue(self, res):
                return res.id[0] == " " and res.get_resname() in _STD_AA

        io = PDBIO()
        io.set_structure(st)
        io.save(out, ProtSelect())
        return out

    # PDB: filtrado por líneas, rápido y sin perder formato
    with open(path) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                if resname in _STD_AA:
                    fout.write(line)
            elif line.startswith("TER"):
                fout.write(line)
        fout.write("END\n")
    return out


# ───────────────────────── PyRosetta: ΔG ──────────────────────────────────────
def compute_rosetta_dg(model_pdb, relax_mode="restrained"):
    """ΔG (REU) con ref2015. relax_mode: 'none' | 'restrained' | 'full'."""
    import pyrosetta
    from pyrosetta import pose_from_file, get_fa_scorefxn
    from pyrosetta.rosetta.protocols.relax import FastRelax

    pyrosetta.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false")

    # Limpieza: el ΔG debe medir SOLO la cadena proteica, igual para el original
    # (que puede traer Zn, ligandos, glicerol, aguas...) que para los modelos
    # (que son proteína sola). Si no se limpia, el ΔΔG queda contaminado por los
    # heteroátomos. Nos quedamos solo con líneas ATOM (descarta HETATM/HOH/ZN/etc.).
    clean_pdb = _strip_hetatm(model_pdb)

    pose = pose_from_file(clean_pdb)
    sfxn = get_fa_scorefxn()  # ref2015 por defecto

    if relax_mode == "none":
        return float(sfxn(pose))

    fr = FastRelax()
    fr.set_scorefxn(sfxn)

    if relax_mode == "restrained":
        # Restricciones de coordenadas sobre todos los átomos (ancla la geometría)
        from pyrosetta.rosetta.protocols.constraint_generator import (
            CoordinateConstraintGenerator, AddConstraints
        )
        cg = CoordinateConstraintGenerator()
        ac = AddConstraints()
        ac.add_generator(cg)
        ac.apply(pose)
        # activa el término de constraint en la score function
        from pyrosetta.rosetta.core.scoring import coordinate_constraint
        sfxn.set_weight(coordinate_constraint, 1.0)
        fr.constrain_relax_to_start_coords(True)

    fr.apply(pose)
    return float(sfxn(pose))


# ───────────────────────── pLDDT (B-factor) ───────────────────────────────────
def compute_mean_plddt(model_pdb):
    """Media del B-factor de los CA (pLDDT por residuo). Lee el PDB con biopython."""
    parser = get_parser(model_pdb)
    s = parser.get_structure("m", model_pdb)
    model = next(s.get_models())
    vals = []
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            if "CA" in res:
                vals.append(res["CA"].get_bfactor())
    return float(np.mean(vals)) if vals else float("nan")


# ───────────────────────── nº cisteínas ───────────────────────────────────────
def compute_n_cysteines(model_pdb):
    parser = get_parser(model_pdb)
    s = parser.get_structure("m", model_pdb)
    model = next(s.get_models())
    n = 0
    for chain in model:
        for res in chain:
            if res.id[0] == " " and res.get_resname() == "CYS":
                n += 1
    return n


# ───────────────────────── carga neta (PROPKA) ────────────────────────────────
def compute_net_charge(model_pdb, ph=7.0):
    """Carga neta a pH dado usando pKa por entorno (PROPKA)."""
    try:
        import propka.run as pk
        mol = pk.single(model_pdb, optargs=["--quiet"])
        # propka expone la carga vs pH; tomamos el valor al pH pedido
        # getPI / charge profile: usamos la suma sobre grupos ionizables
        groups = mol.conformations['AVR'].groups
        charge = 0.0
        for g in groups:
            if not g.titratable:
                continue
            pka = g.pka_value
            q = g.charge  # +1 (base) o -1 (ácido) en estado protonado/desprotonado
            # fracción protonada (Henderson-Hasselbalch)
            if q > 0:   # grupo básico: carga + cuando pH < pKa
                frac = 1.0 / (1.0 + 10 ** (ph - pka))
                charge += frac
            else:       # grupo ácido: carga - cuando pH > pKa
                frac = 1.0 / (1.0 + 10 ** (pka - ph))
                charge -= frac
        return float(charge)
    except Exception as e:
        sys.stderr.write(f"[WARN] PROPKA falló: {e}\n")
        return float("nan")


# ───────────────────────── SAP (solubilidad) ──────────────────────────────────
# Spatial Aggregation Propensity (Chennamsetty et al. 2009), versión simplificada:
#   por cada átomo de cadena lateral, suma de (SASA_relativa * hidrofobicidad)
#   de los residuos cuyo centro cae dentro de un radio R; promedio sobre la proteína.
# Aproximación a nivel residuo con SASA por freesasa + escala de hidrofobicidad.
KD = {  # Kyte-Doolittle (hidrofobicidad; + = hidrofóbico)
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}
# SASA máxima por residuo (tripéptido Gly-X-Gly, Tien 2013 theoretical)
MAXASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

def compute_sap(model_pdb, radius=5.0):
    """SAP score promedio. Mayor = más propenso a agregar = peor solubilidad."""
    import freesasa

    # SASA por residuo con freesasa (requiere PDB; si es CIF, convertir antes)
    struct = freesasa.Structure(model_pdb)
    result = freesasa.calc(struct)
    res_areas = result.residueAreas()  # dict chain -> resnum(str) -> ResidueArea

    parser = get_parser(model_pdb)
    s = parser.get_structure("m", model_pdb)
    model = next(s.get_models())

    # centro de cada residuo (centroide de la cadena lateral, o CA si no hay)
    residues = []
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            rn = res.get_resname()
            if rn not in KD:
                continue
            sc = [a.get_coord() for a in res if a.get_name() not in ("N", "C", "O", "CA")]
            center = np.mean(sc, axis=0) if sc else (
                res["CA"].get_coord() if "CA" in res else None)
            if center is None:
                continue
            # SASA relativa
            try:
                area = res_areas[chain.id][str(res.id[1])].total
            except (KeyError, AttributeError):
                area = 0.0
            rel_sasa = min(area / MAXASA[rn], 1.0) if MAXASA[rn] > 0 else 0.0
            residues.append({
                "center": np.array(center),
                "rel_sasa": rel_sasa,
                "hydro": KD[rn],
            })

    if not residues:
        return float("nan")

    centers = np.array([r["center"] for r in residues])
    sap_per_res = []
    for i, ri in enumerate(residues):
        d = np.linalg.norm(centers - ri["center"], axis=1)
        neigh = d <= radius
        # SAP del residuo i = suma sobre vecinos de (rel_sasa * hydrofobicidad)
        sap_i = sum(residues[j]["rel_sasa"] * residues[j]["hydro"]
                    for j in np.where(neigh)[0])
        sap_per_res.append(sap_i)

    return float(np.mean(sap_per_res))


# ───────────────────────── TM-score ───────────────────────────────────────────
# Mapa de 3 letras -> 1 letra (tmtools necesita la secuencia)
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _ca_coords_and_seq(path, parser):
    """Devuelve (coords Nx3 de CA, secuencia 1-letra) del primer modelo."""
    s = parser.get_structure("s", path)
    model = next(s.get_models())
    coords, seq = [], []
    for chain in model:
        for res in chain:
            if res.id[0] != " " or "CA" not in res:
                continue
            rn = res.get_resname()
            if rn not in THREE_TO_ONE:
                continue
            coords.append(res["CA"].get_coord())
            seq.append(THREE_TO_ONE[rn])
    return np.array(coords, dtype=float), "".join(seq)


def compute_tm_score(model_pdb, ref_pdb):
    """TM-score del modelo vs referencia con tmtools (algoritmo de Zhang:
    rotación óptima que maximiza el TM-score, con alineamiento por secuencia).
    Devuelve el TM normalizado por la longitud de la REFERENCIA (tm_norm_chain2).
    """
    from tmtools import tm_align

    # parser propio para cada fichero (modelo y ref pueden tener extensión distinta)
    coords_m, seq_m = _ca_coords_and_seq(model_pdb, get_parser(model_pdb))
    coords_r, seq_r = _ca_coords_and_seq(ref_pdb,   get_parser(ref_pdb))

    if len(coords_m) < 3 or len(coords_r) < 3:
        sys.stderr.write("[WARN] <3 CA; TM-score=NA\n")
        return float("nan")

    res = tm_align(coords_m, coords_r, seq_m, seq_r)
    # tm_norm_chain2 = normalizado por la longitud de la referencia (chain2 = ref).
    # Es lo coherente con "cuánto recupera el modelo la forma del original".
    return float(res.tm_norm_chain2)


# ───────────────────────── main ───────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--ref", default=None)
    ap.add_argument("--ref-dg", dest="ref_dg", type=float, default=None,
                    help="ΔG del PDB original ya relajado (mismo protocolo). "
                         "Si se da, se calcula rosetta_ddg = rosetta_dg - ref_dg.")
    ap.add_argument("--ref-sap", dest="ref_sap", type=float, default=None,
                    help="SAP del PDB original. Si se da, delta_sap = sap - ref_sap.")
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--relax", default="restrained",
                    choices=["none", "restrained", "full"])
    ap.add_argument("--dg-only", dest="dg_only", action="store_true",
                    help="Solo calcula rosetta_dg y lo escribe (para el PDB original). "
                         "Emite un fichero de una línea con el valor ΔG.")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    # Modo "solo ΔG" (para relajar el original una vez y obtener su ΔG de referencia)
    if args.dg_only:
        # Modo referencia: calcula ΔG Y SAP del original, emite TSV de 2 líneas.
        # (mantiene el flag --dg-only por compatibilidad de nombre)
        out = args.out or f"{args.id}_ref.tsv"
        try:
            dg = compute_rosetta_dg(args.model, args.relax)
        except Exception as e:
            sys.stderr.write(f"[WARN] ΔG original falló: {e}\n")
            dg = float("nan")
        # SAP del original (sobre el PDB tal cual; compute_sap ignora HETATM via KD)
        try:
            sap = compute_sap(args.model)
        except Exception as e:
            sys.stderr.write(f"[WARN] SAP original falló: {e}\n")
            sap = float("nan")
        with open(out, "w") as fh:
            fh.write("metric\tvalue\n")
            fh.write("ref_dg\t" + ("NA" if dg != dg else f"{dg:.6f}") + "\n")
            fh.write("ref_sap\t" + ("NA" if sap != sap else f"{sap:.6f}") + "\n")
        sys.stderr.write(f"[OK] referencia: ΔG={dg}  SAP={sap}\n")
        return

    out = args.out or f"{args.id}_metrics.tsv"
    metrics = {}

    # Si el modelo es CIF (Chai), convertir a PDB para que freesasa/PyRosetta
    # y el resto trabajen homogéneamente.
    model_path = args.model
    if str(model_path).lower().endswith(".cif"):
        try:
            from Bio.PDB import MMCIFParser, PDBIO
            cifp = MMCIFParser(QUIET=True)
            st = cifp.get_structure("m", model_path)
            io = PDBIO()
            io.set_structure(st)
            converted = f"{args.id}_converted.pdb"
            io.save(converted)
            model_path = converted
            sys.stderr.write(f"[INFO] CIF convertido a PDB: {converted}\n")
        except Exception as e:
            sys.stderr.write(f"[WARN] no pude convertir CIF a PDB: {e}\n")

    def safe(name, fn, *a, **k):
        try:
            metrics[name] = fn(*a, **k)
        except Exception as e:
            sys.stderr.write(f"[WARN] {name} falló: {e}\n")
            metrics[name] = float("nan")

    safe("rosetta_dg",  compute_rosetta_dg, model_path, args.relax)
    safe("mean_plddt",  compute_mean_plddt, model_path)
    safe("n_cysteines", compute_n_cysteines, model_path)
    safe("net_charge",  compute_net_charge, model_path, args.ph)
    safe("sap_score",   compute_sap, model_path)
    if args.ref and os.path.exists(args.ref):
        safe("tm_score", compute_tm_score, model_path, args.ref)
    else:
        metrics["tm_score"] = float("nan")

    if "net_charge" in metrics and metrics["net_charge"] == metrics["net_charge"]:
        metrics["net_charge_abs"] = abs(metrics["net_charge"])
    else:
        metrics["net_charge_abs"] = float("nan")

    # ΔΔG vs original: rosetta_dg(modelo) - ref_dg(original).
    # Más negativo = el diseño es MÁS estable que el original.
    dg = metrics.get("rosetta_dg")
    if args.ref_dg is not None and dg is not None and dg == dg:  # dg==dg descarta NaN
        metrics["rosetta_ddg"] = float(dg) - float(args.ref_dg)
    else:
        metrics["rosetta_ddg"] = float("nan")

    # ΔSAP vs original: sap(modelo) - ref_sap(original).
    # Más negativo = el diseño es MENOS agregante = MÁS soluble que el original.
    sap = metrics.get("sap_score")
    if args.ref_sap is not None and sap is not None and sap == sap:
        metrics["delta_sap"] = float(sap) - float(args.ref_sap)
    else:
        metrics["delta_sap"] = float("nan")

    with open(out, "w") as fh:
        fh.write("id\tmetric\tvalue\n")
        for k, v in metrics.items():
            vs = "NA" if (v != v) else f"{v:.6g}"  # v!=v detecta NaN
            fh.write(f"{args.id}\t{k}\t{vs}\n")

    sys.stderr.write(f"[OK] {args.id}: " +
                     ", ".join(f"{k}={metrics[k]}" for k in metrics) + "\n")


if __name__ == "__main__":
    main()
