
# ------------------------------------------------------------
# Generate selective minimization command from a .geo file
# ------------------------------------------------------------
# Usage:
#   ./generate_minimize_bond_lengths_sel.sh model.geo
# ------------------------------------------------------------

if [ $# -ne 1 ]; then
  echo "Usage: $0 model.geo"
  exit 1
fi

GEO="$1"
if [ ! -f "$GEO" ]; then
  echo "Error: file '$GEO' not found."
  exit 1
fi

# Derive associated model name safely
if [[ "$GEO" == *.pdb.geo ]]; then
  MODEL="${GEO%.pdb.geo}.pdb"
else
  MODEL="${GEO%.geo}.pdb"
fi

echo "Processing geometry file: $GEO"
echo "Associated model:         $MODEL"
echo

# 1 Extract |Z| ≥ 4 bond-length outliers
awk '
function abs(x){return x<0?-x:x}
# start of a bond block
/^bond pdb=/ {
  a1=$0; getline; a2=$0;          # two atom lines
  getline;                         # skip "ideal model..." header
  getline; num=$0;                 # numeric values line
  gsub(/^ +| +$/,"",num)
  split(num,f,/ +/)
  ideal=f[1]+0; model=f[2]+0; sigma=f[4]+0
  if (sigma>0) {
    z=abs((model-ideal)/sigma)
    if (z>=4) {
      # extract only atom selectors
      sub(/.*pdb="/,"",a1); sub(/".*/,"",a1)
      sub(/.*pdb="/,"",a2); sub(/".*/,"",a2)
      printf "Z=%6.2f  %s | %s  ideal=%s  model=%s  sigma=%s\n",
             z, a1, a2, f[1], f[2], f[4]
    }
  }
}' "$GEO" | sort -r -k1,1 > bond_outliers.txt

NOUT=$(wc -l < bond_outliers.txt)
echo "Found $NOUT bond length outliers."
if [ "$NOUT" -eq 0 ]; then
  echo "No significant outliers. Exiting."
  exit 0
fi

# 2 Extract unique chain+resid pairs
awk -F'\\|' '
{
  for (i=1;i<=2;i++) {
    s=$i
    gsub(/^[ \t]+|[ \t]+$/, "", s)
    if (match(s, /[ \t]+[A-Za-z0-9*]+[ \t]+([^ \t])[ \t]+([0-9]+)/, m)) {
      key=m[1]" "m[2]; seen[key]=1
    }
  }
}
END {
  for (k in seen) print k
}' bond_outliers.txt \
| sort -k1,1 -k2,2n \
| awk '{printf "chain %s and resid %s\n", $1, $2}' \
> minimize_selection.txt

NRES=$(wc -l < minimize_selection.txt)
echo "Extracted $NRES unique residues to minimize."
echo

# 3 Build Phenix selection string
SEL="$(awk '{printf "(%s) or ", $0}' minimize_selection.txt | sed 's/ or $//')"

# 4 Echo ready-to-run command
echo "✅ Suggested minimization command:"
echo
echo "phenix.geometry_minimization $MODEL selection=\"$SEL\" max_iterations=500"
echo
echo "Files written:"
echo "  - bond_outliers.txt"
echo "  - minimize_selection.txt"
