# @Author: guanja
# @Date:   2019-08-06 20:49:46
# @Last Modified by:   guanja
# @Last Modified time: 2019-08-06 20:51:34


ROOT="/home/guanja/projects/cpkgs/TaroneCpp/TaroneCpp"

EXEC="${ROOT}/compiled/edge_epistasis"

DAT_FN="${ROOT}/examples/data_ihgc/ma_mo_maf_0.05_X.txt"
COV_FN="${ROOT}/examples/data_ihgc/ma_mo.hg19.union_ncov_3.txt"
EDG_FN="${ROOT}/examples/data_ihgc/InWeb3_HUGO_symmetric.tsv"
LAB_FN="${ROOT}/examples/data_ihgc/ma_mo_maf_0.05_Y.txt"
MAP_FN="${ROOT}/examples/data_ihgc/ma_mo_maf_0.05_snp_map.txt"
SNP_FN="${ROOT}/examples/data_ihgc/ma_mo_maf_0.05_snpID.txt"

N_THREADS=2

OUT_PR="${ROOT}/examples/output_ihgc/edgeEpi_ma_mo"

"${EXEC}" \
  -i "${DAT_FN}" \
  -l "${LAB_FN}" \
  -c "${COV_FN}" \
  -m "${MAP_FN}" \
  -e "${EDG_FN}" \
  -s "${SNP_FN}" \
  -f 0.05 \
  -r \
  -n "${N_THREADS}" \
  -o "${OUT_PR}" \