# Selective Geometry Minimization 🧩  
### Targeted bond-angle cleanup in Phenix  

This script semi-automates identification and repair of **bond-angle outliers** in refined cryo-EM or crystallographic RNA models — minimizing only the problematic residues rather than relaxing the entire structure.  

It reads the `.geo` restraint report from Phenix, extracts residues with large |Z| scores, and generates a ready-to-run **phenix.geometry_minimization** command targeting those sites.  

---

## Requirements ⚙️  

- **Bash** (Linux/macOS shell)
- **Python 3.7+** 
- **Phenix ≥ 1.20** (for geometry minimization & restraint generation)  
- Standard Unix tools: `awk`, `sed`, `sort`, `wc`  

---

## Purpose 🧠  

Identify and selectively minimize **bond-angle outliers** (|Z| ≥ 4) from a Phenix `.geo` file, improving local geometry without distorting well-refined regions.  

---

## Protocol 🚀  

### **1. Generate the geometry restraints files (.geo & .eff)**  

From your refined model:  

```bash
phenix.pdb_interpretation yourmodel.pdb write_geo_file=True
```

From a fasta file and dot-bracket secondary structure:
```bash
python dot_bracket_to_eff.py fasta.fasta secstruct.txt \
    --chain A --format phenix > SS.eff
```

### **2. Extract residues with bad angles and generate the geometry-minimization command**

Run the script:  

```bash
generate_minimize_sel.sh yourmodel.pdb.geo
```
This will:
- Parse the `.geo` file  
- Identify all angle outliers (|Z| ≥ 4)  
- Collect the corresponding chain and residue numbers  

Output:
- `angle_outliers.txt` — full list of outliers  
- `minimize_selection.txt` — residue-selection lines  


Print the suggested command:
```
phenix.geometry_minimization yourmodel.pdb \
    selection="(chain A and resid 153) or (chain A and resid 154)" \
    max_iterations=500
```
Optionally add the .eff file to the command to enforce base WC base pairs

### **3. Acknowledgment** 📜

If this script saves you time or sanity, buy me a coffee sometime ☕

**McRae, E.K.S. (2025)** , 
Center for RNA Therapeutics, Houston Methodist Research Institute
