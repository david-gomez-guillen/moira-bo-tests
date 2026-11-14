"""The Python half of MOIRA-BO: the Gaussian processes of the outputs and the
inverse model. moira_bo.R loads it through reticulate and does the rest.

A Gaussian process is fitted per output of the model, on parameters scaled to
[0, 1]. Its kernel is the anisotropic squared exponential, with no noise term:
the models are deterministic or seeded, and a tiny nugget only keeps the
Cholesky factor stable. Its hyperparameters are fitted to the marginal
likelihood by L-BFGS-B, from several starts in a full fit and from the values
it already has in a warm one.

The inverse model is a conditional variational autoencoder of the parameters
given the misses of the outputs. Sampled with every miss at 0, it proposes
parameter vectors that meet every target.
"""

import math

import numpy as np
import torch
from scipy.optimize import minimize
from torch import nn

torch.set_default_dtype(torch.float64)


def set_seed(seed):
    """Seeds numpy and torch, which R's set.seed() does not reach."""
    np.random.seed(int(seed) % (2**32))
    torch.manual_seed(int(seed))


class GP:
    """A Gaussian process with an anisotropic squared exponential kernel,
    k(x, x') = s^2 exp(-sum_i (x_i - x'_i)^2 / (2 l_i^2)), on standardized y.
    """

    NUGGET = 1e-9
    LOG_LS_BOUNDS = (-4.0, 2.0)
    LOG_OS_BOUNDS = (-3.0, 3.0)

    def __init__(self, d):
        self.d = int(d)
        self.theta = np.zeros(self.d + 1)  # log lengthscales, then log scale
        self.constant = True
        self.y_mean = 0.0
        self.y_std = 1.0

    def _kernel(self, X1, X2, theta):
        ls2 = torch.exp(2.0 * theta[: self.d])
        os2 = torch.exp(2.0 * theta[self.d])
        sq = (((X1.unsqueeze(1) - X2.unsqueeze(0)) ** 2) / ls2).sum(-1)
        return os2 * torch.exp(-0.5 * sq)

    def _chol(self, X, theta):
        K = self._kernel(X, X, theta)
        eye = torch.eye(X.shape[0])
        for jitter in [0.0] + [10.0**e for e in range(-8, 1)]:
            L, info = torch.linalg.cholesky_ex(K + (self.NUGGET + jitter) * eye)
            if int(info) == 0:
                return L
        raise RuntimeError("the kernel matrix is not positive definite")

    def _nll(self, theta_np, X, y):
        theta = torch.tensor(theta_np, requires_grad=True)
        try:
            L = self._chol(X, theta)
        except RuntimeError:
            return math.inf, np.zeros_like(theta_np)
        alpha = torch.cholesky_solve(y.unsqueeze(1), L).squeeze(1)
        loss = 0.5 * (y @ alpha) + torch.log(torch.diagonal(L)).sum()
        if not torch.isfinite(loss):
            return math.inf, np.zeros_like(theta_np)
        loss.backward()
        return float(loss.detach()), theta.grad.numpy().copy()

    def fit(self, X, y, starts=5, n_iter=100, warm=False):
        """Fits the hyperparameters: from `starts` points, the first the values
        held, or from the values held alone when `warm`."""
        X = np.asarray(X, dtype=np.float64).reshape(-1, self.d)
        y = np.asarray(y, dtype=np.float64).ravel()
        self.X = torch.tensor(X)
        if y.size < 2 or np.std(y) < 1e-12:
            self.constant = True
            self.y_mean = float(np.mean(y)) if y.size else 0.0
            return self
        self.y_mean = float(np.mean(y))
        self.y_std = float(max(np.std(y), 1e-10))
        Y = torch.tensor((y - self.y_mean) / self.y_std)
        bounds = [self.LOG_LS_BOUNDS] * self.d + [self.LOG_OS_BOUNDS]
        x0s = [self.theta.copy()]
        if not (warm and not self.constant):
            x0s += [np.concatenate([np.random.uniform(-2.0, 1.0, self.d),
                                    np.random.uniform(-1.0, 1.0, 1)])
                    for _ in range(int(starts) - 1)]
        best, best_loss = self.theta.copy(), math.inf
        for x0 in x0s:
            try:
                res = minimize(self._nll, x0, args=(self.X, Y), jac=True, method="L-BFGS-B",
                               bounds=bounds,
                               options={"maxiter": int(n_iter), "ftol": 1e-12, "gtol": 1e-8})
                if np.isfinite(res.fun) and res.fun < best_loss:
                    best, best_loss = res.x.copy(), float(res.fun)
            except Exception:
                pass
        self.theta = best
        theta = torch.tensor(self.theta)
        try:
            L = self._chol(self.X, theta)
        except RuntimeError:
            self.constant = True
            return self
        self.L = L
        self.alpha = torch.cholesky_solve(Y.unsqueeze(1), L).squeeze(1)
        self.constant = False
        return self

    def predict(self, X_new):
        """Mean and standard deviation at the rows of X_new, in y's units."""
        X_new = torch.tensor(np.asarray(X_new, dtype=np.float64).reshape(-1, self.d))
        if self.constant:
            n = X_new.shape[0]
            return np.full(n, self.y_mean), np.zeros(n)
        theta = torch.tensor(self.theta)
        with torch.no_grad():
            Ks = self._kernel(X_new, self.X, theta)
            mu = Ks @ self.alpha
            v = torch.linalg.solve_triangular(self.L, Ks.T, upper=False)
            var = torch.exp(2.0 * theta[self.d]) - (v * v).sum(0)
            sd = torch.sqrt(torch.clamp(var, min=0.0))
        return (mu * self.y_std + self.y_mean).numpy(), (sd * self.y_std).numpy()

    def log_lengthscales(self):
        return self.theta[: self.d].copy()


def fit_gps(gps, X, Y, full, starts=5, n_iter=100):
    """One GP per column of Y, fitted afresh when `full` or when there are none
    yet, and from the hyperparameters of `gps` otherwise. Returns the list."""
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64).reshape(X.shape[0], -1)
    if gps is None or full or len(gps) != Y.shape[1]:
        gps = [GP(X.shape[1]) for _ in range(Y.shape[1])]
    return [gp.fit(X, Y[:, k], starts=starts, n_iter=n_iter, warm=not full)
            for k, gp in enumerate(gps)]


def predict_gps(gps, X):
    """The means and standard deviations of every GP at the rows of X, as two
    matrices of a column per output."""
    preds = [gp.predict(X) for gp in gps]
    return (np.column_stack([p[0] for p in preds]), np.column_stack([p[1] for p in preds]))


def mean_log_lengthscales(gps):
    return np.mean([gp.log_lengthscales() for gp in gps], axis=0)


class _CVAE(nn.Module):
    def __init__(self, d_x, d_c, d_z, hidden):
        super().__init__()
        self.d_z = d_z
        self.enc = nn.Sequential(nn.Linear(d_x + d_c, hidden), nn.ReLU(),
                                 nn.Linear(hidden, hidden), nn.ReLU())
        self.enc_mu = nn.Linear(hidden, d_z)
        self.enc_lv = nn.Linear(hidden, d_z)
        self.dec = nn.Sequential(nn.Linear(d_z + d_c, hidden), nn.ReLU(),
                                 nn.Linear(hidden, hidden), nn.ReLU(),
                                 nn.Linear(hidden, d_x))

    def forward(self, x, c):
        h = self.enc(torch.cat([x, c], 1))
        mu, lv = self.enc_mu(h), self.enc_lv(h).clamp(-10, 5)
        z = mu + torch.randn_like(mu) * torch.exp(0.5 * lv)
        return self.dec(torch.cat([z, c], 1)), mu, lv


class InverseModel:
    """A CVAE with the standardization of the data it was first trained on,
    which later rounds keep, so that a condition means the same across them."""

    def __init__(self, X, C, latent_dim, hidden):
        self.net = _CVAE(X.shape[1], C.shape[1], int(latent_dim), int(hidden))
        self.x_mean, self.x_sd = X.mean(0), X.std(0) + 1e-12
        self.c_mean, self.c_sd = C.mean(0), C.std(0) + 1e-12


def train_cvae(X, C, W, latent_dim=4, hidden=64, epochs=300, lr=1e-3, beta=0.05,
               batch_size=128, seed=0, model=None):
    """Trains an inverse model of the rows of X given those of C, each row's
    loss weighted by W, or keeps training `model`. Returns the model."""
    torch.manual_seed(int(seed))
    X = np.asarray(X, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64).reshape(X.shape[0], -1)
    W = np.asarray(W, dtype=np.float64).ravel()
    if model is None:
        model = InverseModel(X, C, latent_dim, hidden)
    net = model.net
    xs = torch.tensor((X - model.x_mean) / model.x_sd)
    cs = torch.tensor((C - model.c_mean) / model.c_sd)
    ws = torch.tensor(W / W.mean())
    opt = torch.optim.Adam(net.parameters(), lr=float(lr))
    n = X.shape[0]
    net.train()
    for _ in range(int(epochs)):
        perm = torch.randperm(n)
        for i in range(0, n, int(batch_size)):
            idx = perm[i:i + int(batch_size)]
            xr, mu, lv = net(xs[idx], cs[idx])
            rec = ((xr - xs[idx]) ** 2).sum(1)
            kl = -0.5 * (1 + lv - mu**2 - lv.exp()).sum(1)
            loss = (ws[idx] * (rec + float(beta) * kl)).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
    net.eval()
    return model


def sample_cvae(model, c, n, seed=0):
    """n parameter vectors for the condition c, from the prior of the latent."""
    torch.manual_seed(int(seed))
    net = model.net
    c = (np.asarray(c, dtype=np.float64).ravel() - model.c_mean) / model.c_sd
    with torch.no_grad():
        z = torch.randn(int(n), net.d_z)
        x = net.dec(torch.cat([z, torch.tensor(np.tile(c, (int(n), 1)))], 1)).numpy()
    return x * model.x_sd + model.x_mean
