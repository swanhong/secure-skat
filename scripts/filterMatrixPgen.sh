#!/bin/sh
PGEN_PREFIX=$1
NR=$2
NC=$3
ROW_FILT_FILE=$4
COL_FILT_FILE=$5
COL_NAMES_FILE=$6
COL_START_POS=$7
OUT_FILE=$8
SECURE_ORIENTATION=${9:-false}

python3 scripts/filterLines.py "${COL_NAMES_FILE}" "${COL_FILT_FILE}" "${COL_START_POS}" "${OUT_FILE}.snps.txt"
if [ -z "${ROW_FILT_FILE}" ]
then
    plink2 --pfile "${PGEN_PREFIX}" --extract "${OUT_FILE}.snps.txt" --indiv-sort none --make-bed --out "${OUT_FILE}"
else
    plink2 --pfile "${PGEN_PREFIX}" --keep "${ROW_FILT_FILE}" --extract "${OUT_FILE}.snps.txt" --indiv-sort none --make-bed --out "${OUT_FILE}"
fi
EXTRA_ARGS=""
if [ "${SECURE_ORIENTATION}" = "true" ]
then
    EXTRA_ARGS="--secure-orientation"
fi
python3 scripts/plinkBedToBinary.py "${OUT_FILE}.bed" "${NR}" "${NC}" "${OUT_FILE}" ${EXTRA_ARGS}
