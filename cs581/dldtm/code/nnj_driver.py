"""Run NeuralNJ (Zhang et al., MBE 2025, doi:10.1093/molbev/msaf260) argmax inference
(no RL search, no raxml branch-length optimisation) on PHYLIP files, CPU only.
Runs inside the Python 3.10 env (/opt/mm/root/envs/nnj) with the pure-Python raxmlpy shim.
Usage: python nnj_driver.py IN1.phy OUT1.tre [IN2.phy OUT2.tre ...]
"""
import os
import sys

NNJ = os.environ.get("NNJ_REPO", "/opt/src/NeuralNJ")
sys.path.insert(0, NNJ)
sys.path.insert(0, os.environ.get("RAXMLPY_SHIM", "/opt/src/raxmlpy_shim"))
os.chdir(NNJ)
import torch  # noqa: E402
import utils  # noqa: E402
import finetune_rl_search as F  # noqa: E402

torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
cfgs = utils.empty_config()
cfgs.merge_from_file(f"{NNJ}/config/finetune_reinforce_search_example.yaml")
cfgs.env.batch_size = 1
utils.set_evolution_model("GTR+I+G")
F.cfgs = cfgs
F.evolution_model = "GTR+I+G"
env = F.PhyInferEnv(cfgs, F.device)
model = F.PGPI(cfgs).to(F.device)
ck = torch.load(f"{NNJ}/checkpoint/train_on_generateemprical_GTR+I+G_final.pt", map_location="cpu", weights_only=False)
model.load_state_dict(ck["model_state_dict"])
model.eval()
args = sys.argv[1:]
for phy, out in zip(args[::2], args[1::2]):
    with torch.no_grad():
        r = F.Agmax_one_instance(cfgs, phy, model, env, branch_optimize=False)
    with open(out, "w") as f:
        f.write(r["best_tree_str"].strip() + "\n")
