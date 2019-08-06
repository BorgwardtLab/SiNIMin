# @Author: guanja
# @Date:   2019-07-10 14:50:08
# @Last Modified by:   guanja
# @Last Modified time: 2019-07-31 20:03:27


ROOT="/home/guanja/projects/cpkgs/TaroneCpp/TaroneCpp"

EXEC="${ROOT}/compiled/edge_epistasis"

DAT_FN="${ROOT}/examples/data/data.txt"
COV_FN="${ROOT}/examples/data/covar.txt"
EDG_FN="${ROOT}/examples/data/edges.txt"
LAB_FN="${ROOT}/examples/data/labels.txt"
MAP_FN="${ROOT}/examples/data/mapping_file.txt"
SNP_FN="${ROOT}/examples/data/snp_ids.txt"

N_THREADS=2

OUT_PR="${ROOT}/examples/output/edgeEpiCMH_examples"

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

