"""Re-predict multi-label models on chr 1-7 and save averaged LOGITS (forward + strand-flipped), not sigmoid probabilities,
so plateaus where p saturates at 1.0 keep their internal structure. Output exp_ml/logit_<name>.npz, keys 'pool__<chrom>', (2,3,n20) float16."""
import sys, os, numpy as np, torch
sys.path[:0] = ['.', 'boundary_model']
os.environ.setdefault('POOL', '1')
import pfprep as P, seg_features as SF, exp as X_, exp_ml as XM
name = sys.argv[1]; torch.set_num_threads(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
ck = torch.load(f'exp_ml/{name}.pt', weights_only=False); a = ck['args']
spec = set(a['spec'].split(',')); TRK = [x for x in a['tracks'].split(',') if x] or None
X20, nu, kg, oh = SF.assemble(spec - {'strand', 'daft'}, TRK)
if spec & {'strand', 'daft'}:
    ST = X_.stranded(spec); X_.NSTR = next(iter(ST.values())).shape[0] // 2
    for k in X20: n = min(X20[k].shape[1], ST[k].shape[1]); X20[k] = np.concatenate([X20[k][:, :n], ST[k][:, :n]])
keys = [k for k in X20 if os.environ.get('PRED_ALL') or int(k.split('|')[1][6:8]) <= 7]
X = {k: X_.to_res(X20[k], a['res']) for k in keys}
net = X_.make_net(a, ck['n_in']); net.load_state_dict(ck['state']); net.eval()
win = a['win_bp'] // a['res']; win -= win % 8; crop = max(8, (a['crop_bp'] // a['res']) // 8 * 8)
out = {}
with torch.no_grad():
    for k in keys:
        x = X[k]; nb = x.shape[1]; acc = np.zeros((6, nb), np.float32); cnt = np.zeros(nb, np.float32)
        step = win - 2 * crop; starts = list(range(0, max(1, nb - win + 1), step))
        if starts[-1] + win < nb: starts.append(max(0, nb - win))
        for s in starts:
            xb = torch.from_numpy(x[:, s:s + win].astype(np.float32))[None]; lg = XM.logits3(net, xb)
            d = torch.zeros(1, 2, 3, xb.shape[-1]); xf, _, _ = XM.flip3(xb, d, d, nu, kg, oh)
            lg = 0.5 * (lg + XM.logits3(net, xf).flip(-1)[:, [1, 0]])
            p = lg.reshape(6, -1).numpy(); L = p.shape[1]; lo = 0 if s == 0 else crop; hi = L if s + win >= nb else L - crop
            acc[:, s + lo:s + hi] += p[:, lo:hi]; cnt[s + lo:s + hi] += 1
        lgt = (acc / np.maximum(cnt, 1)).reshape(2, 3, nb)
        n20 = X20[k].shape[1]; out[k.replace('|', '__')] = X_.pred_to20(lgt, a['res'], n20).astype(np.float16)
np.savez_compressed(f"exp_ml/{'logitall' if os.environ.get('PRED_ALL') else 'logit'}_{name}.npz", **out); print('ok', name, len(out))
