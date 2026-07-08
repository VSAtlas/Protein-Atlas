#!/usr/bin/env bash
set -euo pipefail

MIN_ACT=30
MAX_ACT=50
DECOYS_PER_ACT=40
MAX_DEC=3000

for lib in pur2 mk01 fabp4; do
  src_pdbqt="prepped_ligands/${lib}"
  src_sdf="extracted_ligands/${lib}"

  dst_pdbqt="prepped_ligands/bench_${lib}"
  dst_sdf="extracted_ligands/bench_${lib}"

  mkdir -p "${dst_pdbqt}/actives" "${dst_pdbqt}/decoys" "${dst_sdf}"

  shopt -s nullglob
  act_files=( "${src_pdbqt}"/actives_final_*.pdbqt )
  dec_files=( "${src_pdbqt}"/decoys_final_*.pdbqt )

  n_act=${#act_files[@]}
  n_dec=${#dec_files[@]}

  if [ "$n_act" -eq 0 ] || [ "$n_dec" -eq 0 ]; then
    echo "[bench] ${lib}: missing actives_final_*.pdbqt or decoys_final_*.pdbqt in ${src_pdbqt}"
    continue
  fi

  take_act=$(( (n_act + 9) / 10 ))
  if [ "$take_act" -lt "$MIN_ACT" ] && [ "$n_act" -ge "$MIN_ACT" ]; then take_act="$MIN_ACT"; fi
  if [ "$take_act" -gt "$MAX_ACT" ]; then take_act="$MAX_ACT"; fi
  if [ "$take_act" -gt "$n_act" ]; then take_act="$n_act"; fi

  take_dec=$(( (n_dec + 9) / 10 ))
  min_dec=$(( take_act * DECOYS_PER_ACT ))
  if [ "$take_dec" -lt "$min_dec" ] && [ "$n_dec" -ge "$min_dec" ]; then take_dec="$min_dec"; fi
  if [ "$take_dec" -gt "$MAX_DEC" ]; then take_dec="$MAX_DEC"; fi
  if [ "$take_dec" -gt "$n_dec" ]; then take_dec="$n_dec"; fi

  printf '%s\n' "${act_files[@]}" | head -n "$take_act" | while read -r f; do
    cp --update=none "$f" "${dst_pdbqt}/actives/"
  done

  printf '%s\n' "${dec_files[@]}" | head -n "$take_dec" | while read -r f; do
    cp --update=none "$f" "${dst_pdbqt}/decoys/"
  done

  echo "[bench] ${lib}: PDBQT actives ${take_act}/${n_act}, decoys ${take_dec}/${n_dec} -> ${dst_pdbqt}"

  act_sdf="${src_sdf}/actives_final.sdf"
  dec_sdf="${src_sdf}/decoys_final.sdf"

  if [ -f "$act_sdf" ]; then
    awk -v n="$take_act" 'BEGIN{RS="\\$\\$\\$\\$\\n"; ORS="$$$$\n"} NR<=n {print $0}' \
      "$act_sdf" > "${dst_sdf}/actives_final.sdf"
  else
    echo "[bench] ${lib}: missing ${act_sdf} (skipping SDF actives)"
  fi

  if [ -f "$dec_sdf" ]; then
    awk -v n="$take_dec" 'BEGIN{RS="\\$\\$\\$\\$\\n"; ORS="$$$$\n"} NR<=n {print $0}' \
      "$dec_sdf" > "${dst_sdf}/decoys_final.sdf"
  else
    echo "[bench] ${lib}: missing ${dec_sdf} (skipping SDF decoys)"
  fi

  echo "[bench] ${lib}: SDF actives first ${take_act}, decoys first ${take_dec} -> ${dst_sdf}"
done
