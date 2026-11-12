"""WitnessCell v4 source neural encoder and bounded partner adapter.

SourceEncoderV4 is retrained from fixed 64-D GenePT features; it does not load or
approximate the archived 1536-D checkpoint. Target training uses frozen cached
source embeddings. PartnerHeadV4 is the only target-trainable module. No method
uses test labels, other query predictions, or a query-batch projection.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import torch
from torch import Tensor, nn


@dataclass(frozen=True)
class SourceConfig:
    knowledge_dim: int = 64
    program_dim: int = 32
    width: int = 64
    heads: int = 4
    ff_width: int = 128


@dataclass(frozen=True)
class PartnerConfig:
    width: int = 64
    environment_rank: int = 4
    extra_dim: int = 9
    attention_dim: int = 16
    hidden_dim: int = 16
    bottleneck_dim: int = 8
    dropout: float = 0.1
    use_extra: bool = True


def _expand_vector(x: Tensor, leading: tuple[int, ...], width: int, name: str) -> Tensor:
    if x.shape[-1] != width:
        raise ValueError(f"{name}: expected last dimension {width}, got {tuple(x.shape)}")
    return x.expand(*leading, width)


def _scalar(x: Tensor, leading: tuple[int, ...], name: str) -> Tensor:
    if x.shape == (*leading, 1):
        x = x.squeeze(-1)
    try:
        return x.expand(*leading)
    except RuntimeError as e:
        raise ValueError(f"{name} cannot broadcast to {leading}: {tuple(x.shape)}") from e


class SourceEncoderV4(nn.Module):
    """Two-token RNA-query attention, trained on masked source coordinates.

    pretrain_forward accepts [B,32], [B,64], [B,64], [B,32]. Singleton rows may
    broadcast. The returned [B] value is the source target on the trainer's fixed
    scale. The trainer MUST exclude each predicted coordinate from its input
    response-program score. Mean, PCA and scales must use source training TFs.

    encode_endpoint returns RNA-conditioned z[B,64]; encode_readout returns
    q[B,64]. These two methods support caching before target adaptation.
    """
    def __init__(self, config: SourceConfig | None = None):
        super().__init__()
        self.config = config or SourceConfig()
        c = self.config; d = c.width
        self.endpoint_knowledge = nn.Linear(c.knowledge_dim, d)
        self.endpoint_response = nn.Linear(c.program_dim, d)
        self.token_type = nn.Parameter(torch.randn(2, d) * 0.02)
        self.readout_knowledge = nn.Linear(c.knowledge_dim, d)
        self.readout_loading = nn.Linear(c.program_dim, d, bias=False)
        self.readout_missing = nn.Linear(1, d, bias=False)
        self.readout_norm = nn.LayerNorm(d)
        self.query_norm = nn.LayerNorm(d)
        self.memory_norm = nn.LayerNorm(d)
        self.cross_attention = nn.MultiheadAttention(d, c.heads, dropout=0.0, batch_first=True)
        self.ff_norm = nn.LayerNorm(d)
        self.feedforward = nn.Sequential(nn.Linear(d, c.ff_width), nn.GELU(), nn.Linear(c.ff_width, d))
        self.endpoint_norm = nn.LayerNorm(d)
        self.source_head = nn.Linear(d, 1)

    @staticmethod
    def _rows(*arrays: Tensor):
        if any(a.ndim != 2 for a in arrays):
            raise ValueError("Source features must be two-dimensional")
        n = max(a.shape[0] for a in arrays)
        if any(a.shape[0] not in (1, n) for a in arrays):
            raise ValueError("Source feature row counts must match or be singleton")
        return tuple(a.expand(n, -1) for a in arrays)

    def encode_readout(self, rna_knowledge: Tensor, readout_loadings: Tensor,
                       missing: Tensor | None = None) -> Tensor:
        rna_knowledge, readout_loadings = self._rows(rna_knowledge, readout_loadings)
        n = len(rna_knowledge)
        if missing is None:
            missing = rna_knowledge.new_zeros(n, 1)
        else:
            missing = missing.reshape(-1, 1).expand(n, 1).to(rna_knowledge.dtype)
        return self.readout_norm(self.readout_knowledge(rna_knowledge.detach())
                                 + self.readout_loading(readout_loadings.detach())
                                 + self.readout_missing(missing))

    def encode_endpoint(self, tf_knowledge: Tensor, masked_program: Tensor,
                        readout_query: Tensor) -> Tensor:
        tf_knowledge, masked_program, readout_query = self._rows(tf_knowledge, masked_program, readout_query)
        memory = torch.stack((self.endpoint_knowledge(tf_knowledge.detach()),
                              self.endpoint_response(masked_program)), dim=1)
        memory = self.memory_norm(memory + self.token_type[None, :, :])
        a, _ = self.cross_attention(self.query_norm(readout_query)[:, None, :],
                                    memory, memory, need_weights=False)
        z = readout_query + a[:, 0, :]
        return self.endpoint_norm(z + self.feedforward(self.ff_norm(z)))

    def pretrain_forward(self, masked_program: Tensor, tf_knowledge: Tensor,
                         rna_knowledge: Tensor, readout_loadings: Tensor,
                         missing: Tensor | None = None) -> Tensor:
        q = self.encode_readout(rna_knowledge, readout_loadings, missing)
        z = self.encode_endpoint(tf_knowledge, masked_program, q)
        return self.source_head(z).squeeze(-1)

    def forward(self, masked_program: Tensor, tf_knowledge: Tensor,
                rna_knowledge: Tensor, readout_loadings: Tensor,
                missing: Tensor | None = None) -> Tensor:
        return self.pretrain_forward(masked_program, tf_knowledge, rna_knowledge, readout_loadings, missing)

    def freeze(self):
        self.eval()
        for p in self.parameters():
            p.requires_grad_(False)
            p.grad = None
        return self

    def metadata(self):
        return {"config": asdict(self.config), "parameters": sum(p.numel() for p in self.parameters()),
                "trainable_parameters": sum(p.numel() for p in self.parameters() if p.requires_grad),
                "source_task": "masked-coordinate reconstruction of measured source single-ORF profiles",
                "checkpoint_origin": "fresh 64-D architecture; no 1536-D checkpoint conversion"}


class PartnerHeadV4(nn.Module):
    """An 8,517-parameter partner adapter over frozen neural source features.

    Input shapes share broadcastable leading dimensions, e.g. [E,J,64]:
      z_a,z_b,q: RNA-conditioned frozen endpoint A, endpoint B, and RNA embeddings
      single_a,single_b: centered/scaled single measurements [E,J]
      readout_mean: separately scaled shared mean [J] or [E,J]
      extra_a,extra_b: target-only residual-PC scores [...,9], padded if rank < 9
      extra_readout: corresponding RNA residual-PC loadings [...,9]

    The two original endpoint vectors stay separate throughout the head. A small
    attention query attends jointly to the two adapted endpoint vectors. The
    ordered predictor is evaluated in both endpoint orders, then averaged. This
    enforces exchange symmetry without collapsing endpoint ownership. In eval
    mode, RNA coordinates and query pairs can be evaluated in arbitrary chunks.
    """
    def __init__(self, config: PartnerConfig | None = None):
        super().__init__()
        self.config = config or PartnerConfig()
        c = self.config; d = c.width
        self.environment = nn.Sequential(nn.Linear(d, c.environment_rank), nn.GELU(),
                                         nn.Linear(c.environment_rank, d))
        self.node_extra = nn.Linear(c.extra_dim, d)
        self.readout_extra = nn.Linear(c.extra_dim, d)
        self.attn_q = nn.Linear(d, c.attention_dim)
        self.attn_k = nn.Linear(d, c.attention_dim)
        self.attn_v = nn.Linear(d, c.attention_dim)
        width = 3*d + c.attention_dim + 3
        self.head = nn.Sequential(nn.Linear(width, c.hidden_dim), nn.GELU(), nn.Dropout(c.dropout),
                                  nn.Linear(c.hidden_dim, c.bottleneck_dim), nn.GELU(),
                                  nn.Linear(c.bottleneck_dim, 1))
        for module in (self.environment[-1], self.node_extra, self.readout_extra, self.head[-1]):
            nn.init.zeros_(module.weight)
            nn.init.zeros_(module.bias)

    def _features(self, z_a: Tensor, z_b: Tensor, q: Tensor,
                  single_a: Tensor, single_b: Tensor, readout_mean: Tensor,
                  extra_a: Tensor | None, extra_b: Tensor | None,
                  extra_readout: Tensor | None):
        c = self.config
        leading = torch.broadcast_shapes(z_a.shape[:-1], z_b.shape[:-1], q.shape[:-1])
        za = _expand_vector(z_a.detach(), leading, c.width, "z_a")
        zb = _expand_vector(z_b.detach(), leading, c.width, "z_b")
        query = _expand_vector(q.detach(), leading, c.width, "q")
        za = za + self.environment(za)
        zb = zb + self.environment(zb)
        if c.use_extra:
            if extra_a is not None:
                za = za + self.node_extra(_expand_vector(extra_a.detach(), leading, c.extra_dim, "extra_a"))
            if extra_b is not None:
                zb = zb + self.node_extra(_expand_vector(extra_b.detach(), leading, c.extra_dim, "extra_b"))
            if extra_readout is not None:
                query = query + self.readout_extra(_expand_vector(extra_readout.detach(), leading, c.extra_dim, "extra_readout"))
        qa = self.attn_q(query)
        keys = torch.stack((self.attn_k(za), self.attn_k(zb)), dim=-2)
        values = torch.stack((self.attn_v(za), self.attn_v(zb)), dim=-2)
        weights = ((qa.unsqueeze(-2) * keys).sum(dim=-1) / math.sqrt(c.attention_dim)).softmax(dim=-1)
        attended = (weights.unsqueeze(-1) * values).sum(dim=-2)
        a = _scalar(single_a, leading, "single_a")
        b = _scalar(single_b, leading, "single_b")
        mean = _scalar(readout_mean, leading, "readout_mean")
        f_ab = torch.cat((za, zb, query, attended, torch.stack((a,b,mean),dim=-1)), dim=-1)
        f_ba = torch.cat((zb, za, query, attended, torch.stack((b,a,mean),dim=-1)), dim=-1)
        return f_ab, f_ba, weights

    def forward(self, z_a: Tensor, z_b: Tensor, q: Tensor,
                single_a: Tensor, single_b: Tensor, readout_mean: Tensor,
                extra_a: Tensor | None = None, extra_b: Tensor | None = None,
                extra_readout: Tensor | None = None) -> Tensor:
        f_ab, f_ba, _ = self._features(z_a,z_b,q,single_a,single_b,readout_mean,extra_a,extra_b,extra_readout)
        return 0.5 * (self.head(f_ab) + self.head(f_ba)).squeeze(-1)

    def attention_weights(self, *args, **kwargs) -> Tensor:
        """Diagnostic weights [...,2]; these are not causal explanations."""
        return self._features(*args, **kwargs)[2]

    def metadata(self):
        return {"config": asdict(self.config), "parameters": sum(p.numel() for p in self.parameters()),
                "trainable_parameters": sum(p.numel() for p in self.parameters() if p.requires_grad),
                "initial_prediction": "exactly zero additive correction to unchanged M0",
                "source_features": "detached frozen caches; target loss cannot update source weights",
                "symmetry": "average shared ordered head over (A,B) and (B,A)",
                "target_attention": "one RNA query jointly attends to both complete endpoint vectors",
                "batch_dependent_projection": False,
                "target_extra": "fold-specific residual single-response PCs; scores and RNA loadings padded to 9"}


FrozenSourceEncoder = SourceEncoderV4
PartnerAdapter = PartnerHeadV4

