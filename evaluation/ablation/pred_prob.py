"""Re-create exp_ml/pred_<name>.npz (sigmoid probabilities, chr 6-7 + 1-5) from a saved checkpoint, exactly as exp_ml.py does after training."""
import sys, os, numpy as np, torch
sys.path[:0] = ['.', 'boundary_model']; os.environ.setdefault('POOL', '1')
import seg_features as SF, exp as X_, exp_ml as XM
name = sys.argv[1]; torch.set_num_threads(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
ck = torch.load(f'exp_ml/{name}.pt', weights_only=False); a = ck['args']
spec = set(a['spec'].split(',')); TRK = [x for x in a['tracks'].split(',') if x] or None
X20, nu, kg, oh = SF.assemble(spec - {'strand', 'daft'}, TRK)
if spec & {'strand', 'daft'}:
    ST = X_.stranded(spec); X_.NSTR = next(iter(ST.values())).shape[0] // 2
    for k in X20: n = min(X20[k].shape[1], ST[k].shape[1]); X20[k] = np.concatenate([X20[k][:, :n], ST[k][:, :n]])
Y, W = XM.load_labels(20)
keys = [k for k in X20 if int(k.split('|')[1][6:8]) <= 7]
n20 = {k: min(X20[k].shape[1], Y[k].shape[-1]) for k in keys}
X = {k: X_.to_res(X20[k][:, :n20[k]], a['res']) for k in keys}
net = X_.make_net(a, ck['n_in']); net.load_state_dict(ck['state']); net.eval()
win = a['win_bp'] // a['res']; win -= win % 8; crop = max(8, (a['crop_bp'] // a['res']) // 8 * 8)
Pd = XM.predict(net, X, keys, nu, kg, oh, win, crop); Pd = {k: X_.pred_to20(v, a['res'], n20[k]) for k, v in Pd.items()}
np.savez_compressed(f'exp_ml/pred_{name}.npz', **{k.replace('|', '__'): v for k, v in Pd.items()}); print('ok', name)
