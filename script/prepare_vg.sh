#!/usr/bin/env bash

if [[ $# -ne 1 ]]; then
    echo "[ERROR] Usage: $0 <GRAPH.vg>"
    exit 2
fi

VG_FILE=$(realpath "$1")

echo "[INFO] Preparing '${VG_FILE}'"

echo -n "[INFO] Checking for bgzip... "
BGZIP=$(command -v bgzip)
if [[ -n $BGZIP ]]; then
    echo "found"
else
    echo "no"
fi

echo -n "[INFO] Checking for gzip... "
GZIP=$(command -v gzip)
if [[ -n $GZIP ]]; then
    echo "found"
else
    echo "no"
fi

if [[ -n $BGZIP ]]; then   # Using bgzip
    ZIP_CMD=$BGZIP
    FORMAT="bgzip"
elif [[ -n $GZIP ]]; then  # Using gzip
    ZIP_CMD=$GZIP
    FORMAT="gzip"
else
    echo "[ERROR] No binary found for neither bgzip nor gzip."
    echo
    echo "[HINT] Install bgzip by running one of the following commands:"
    echo "[HINT] →  apt install tabix"
    echo "[HINT] →  conda install tabix"
    echo "[HINT] →  brew install htslib"
    exit 2
fi

if ${ZIP_CMD} -t "${VG_FILE}"; then
    echo "[INFO] No change needed. You're good to go!"
else
    NEW_VG_FILE="${VG_FILE%.vg}_psi_${FORMAT}.vg"
    echo "[INFO] Compressing (${FORMAT}) the input graph..."
    echo "[INFO] New file:  '${NEW_VG_FILE}'"
    if [[ -f ${NEW_VG_FILE} ]]; then
        echo "[ERROR] File exists: ${NEW_VG_FILE}"
        exit 2
    fi
    ${ZIP_CMD} -c "${VG_FILE}" > "${NEW_VG_FILE}"
fi
