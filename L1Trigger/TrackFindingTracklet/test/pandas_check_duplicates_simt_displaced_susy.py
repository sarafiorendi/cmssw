#!/usr/bin/env python3
'''
Check for SIM-level duplicate SimTracks
see issues in https://github.com/cms-sw/cmssw/issues/50253
'''

import sys
import awkward as ak
import numpy as np
import pandas as pd
import uproot
import pdb

MATCH_MODE = "vertex"  # "vertex" (needs st_vx/st_vy/st_vz) or "helix" (d0/z0-based, v1 behaviour)

VTX_TOL = 1e-3   ## cm
DR_TOL = 1e-3    ## angular separation 
DPT_REL_TOL = 1  ## relative pt difference
DZ0_TOL = 1e-3   ## cm 
DD0_TOL = 1e-3   ## cm 

TREE_PATH_CANDIDATES = ["tpNtuple/eventTree", "eventTree"]


def open_tree(path):
    f = uproot.open(path)
    for cand in TREE_PATH_CANDIDATES:
        if cand in f:
            return f[cand]
    raise RuntimeError(f"Could not find eventTree in {path}. Top-level keys: {f.keys()}")


def to_flat_df(tree, branches):
    arrs = tree.arrays(branches, library="ak")
    counts = ak.num(arrs[branches[0]])
    counts_np = ak.to_numpy(counts)
    evt = np.repeat(np.arange(len(counts_np)), counts_np)
    idx = np.concatenate([np.arange(n) for n in counts_np]) if len(counts_np) else np.array([], dtype=int)
    data = {"evt": evt, "idx": idx}
    for b in branches:
        data[b] = ak.to_numpy(ak.flatten(arrs[b]))
    return pd.DataFrame(data), int(len(counts_np))


def find_duplicate_pairs(df):

    counts = df["evt"].value_counts()

    m = df.merge(df, on="evt", suffixes=("_A", "_B"))
    m = m[m["idx_A"] < m["idx_B"]]

    deta = m["eta_A"] - m["eta_B"]
    dphi = (m["phi_A"] - m["phi_B"] + np.pi) % (2 * np.pi) - np.pi
    dR = np.sqrt(deta**2 + dphi**2)
    dpt_rel = (m["pt_A"] - m["pt_B"]).abs() / m["pt_A"].clip(lower=1e-9)
    mask = (
        (dR < DR_TOL)
        & (dpt_rel < DPT_REL_TOL)
#         & ((m["z0_A"] - m["z0_B"]).abs() < DZ0_TOL)
#         & ((m["d0_A"] - m["d0_B"]).abs() < DD0_TOL)
    )
    m = m[mask]
    return m

def main(files):
    total_events = 0
    total_tps = 0
    total_sig_tps = 0
    all_pairs = []
    mode_used = None
    any_match_branch_seen = False

    for fname in files:
        tree = open_tree(fname)
        avail = set(tree.keys())

        base_need = ["st_pt", "st_eta", "st_phi","st_pdgid", "st_parent_pdgid", "st_status", "st_vx", "st_vy", "st_vz"]
        missing_base = [b for b in base_need if b not in avail]
        if missing_base:
            print(
                f"[{fname}] WARNING: missing branches {missing_base}. Available st_* branches: "
                f"{sorted(b for b in avail if b.startswith('st_'))}"
            )
            continue

        need = list(base_need)
        df, nev = to_flat_df(tree, need)
        df = df.rename(columns={b: b.replace("st_", "").replace("matchtrk_", "matchtrk_") for b in need})

        total_events += nev
        total_tps += len(df)

        pairs = find_duplicate_pairs(df)
        pairs = pairs.copy()
        pairs["file"] = fname
        all_pairs.append(pairs)

        print(f"[{fname}] events={nev}  TPs(all)={total_tps}  TPs(signal)={len(df)}  pairs={len(pairs)}")

    pairs = pd.concat(all_pairs, ignore_index=True) if all_pairs else pd.DataFrame()

    print("\n===== SUMMARY =====")
    print(f"Match within dR = {DR_TOL} and relPt = {DPT_REL_TOL}")
    print(f"Total events examined:                  {total_events}")
    print(f"Total SimTracks:                        {total_tps}")
    print(f"Candidate SIM-duplicate SimTrk pairs:   {len(pairs)}")

    if len(pairs):
        mismatched = int((pairs["pdgid_A"] != pairs["pdgid_B"]).sum())
        print(f"  -> with mismatched st_pdgid:            {mismatched} ({100.0 * mismatched / len(pairs):.1f}%)")
        print(f"  duplicate pairs / signal SimTks:        {100.0 * len(pairs) / max(total_sig_tps, 1):.3f}%")

        ## this is when running on L1TF Ntuples
        if any_match_branch_seen and "matchtrk_pt_A" in pairs.columns:
            mi = pairs["matchtrk_pt_A"] > 0
            mj = pairs["matchtrk_pt_B"] > 0
            print(f"  both members matched to an L1 track: {int((mi & mj).sum())}")
            print(f"  exactly one member matched:           {int((mi ^ mj).sum())}  <-- dilutes measured efficiency")
            print(f"  neither member matched:               {int((~mi & ~mj).sum())}")


        print("\n===== DUPLICATE PAIR PROPERTIES (pt, eta, phi, pdgId) =====")
        cols = ["file", "evt", "idx_A", "idx_B", "pt_A", "eta_A", "phi_A", "pdgid_A", "parent_pdgid_A", "status_A", "vx_A", "vy_A", "vz_A", 
                                                 "pt_B", "eta_B", "phi_B", "pdgid_B", "parent_pdgid_B", "status_B", "vx_B", "vy_B", "vz_B"]
        cols = [c for c in cols if c in pairs.columns]
        show = pairs[cols].sort_values(["file", "evt"]).reset_index(drop=True)
        n_show =100
        with pd.option_context("display.max_rows", None, "display.width", 200):
          print(show.head(n_show).to_string(index=False, float_format=lambda x: f"{x:.4f}"))


    else:
        print("No duplicate pairs found. ")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])