# @Author: guanja
# @Date:   2019-07-10 14:50:08
# @Last Modified by:   guanja
# @Last Modified time: 2019-09-17 12:08:59


EXEC="../../SiNIMin/compiled/sinimin"

DAT_FN="../data/example_data_X.txt"
COV_FN="../data/example_data_C.txt"
EDG_FN="../data/example_data_edges.txt"
LAB_FN="../data/example_data_Y.txt"
MAP_FN="../data/example_data_snp_map.txt"
SNP_FN="../data/example_data_snpID.txt"

# Set the parameters for WY-permutations.
N_THREADS=10
N_PERM=1000

OUT_PR="../output"
mkdir -p "${OUT_PR}"

"${EXEC}" \
  -i "${DAT_FN}" \
  -l "${LAB_FN}" \
  -c "${COV_FN}" \
  -m "${MAP_FN}" \
  -e "${EDG_FN}" \
  -s "${SNP_FN}" \
  -f 0.05 \
  -p ${N_PERM} \
  -n ${N_THREADS} \
  -o "${OUT_PR}/example_data_siniminWY" \

