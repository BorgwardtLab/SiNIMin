# @Author: guanja
# @Date:   2019-07-10 14:50:08
# @Last Modified by:   guanja
# @Last Modified time: 2019-09-17 12:10:06


ROOT="/home/guanja/projects/cpkgs/TaroneCpp/TaroneCpp"

EXEC="${ROOT}/compiled/sinimin"

DAT_FN="${ROOT}/examples/simulation/sim_ps_0.05_pcon_0.05_simID_0_X.txt"
COV_FN="${ROOT}/examples/simulation/sim_ps_0.05_pcon_0.05_simID_0_C.txt"
EDG_FN="${ROOT}/examples/simulation/sim_ps_0.05_pcon_0.05_simID_0_edges.txt"
LAB_FN="${ROOT}/examples/simulation/sim_ps_0.05_pcon_0.05_simID_0_Y.txt"
MAP_FN="${ROOT}/examples/simulation/sim_ps_0.05_pcon_0.05_simID_0_snp_map.txt"
SNP_FN="${ROOT}/examples/simulation/sim_ps_0.05_pcon_0.05_simID_0_snpID.txt"

N_THREADS=5

OUT_PR="${ROOT}/examples/output/simulation"
mkdir -p "${OUT_PR}"

"${EXEC}" \
  -i "${DAT_FN}" \
  -l "${LAB_FN}" \
  -c "${COV_FN}" \
  -m "${MAP_FN}" \
  -e "${EDG_FN}" \
  -s "${SNP_FN}" \
  -f 0.05 \
  -d 1 \
  -n 10 \
  -p 1000 \
  -o "${OUT_PR}/edgeEpiWY" \

