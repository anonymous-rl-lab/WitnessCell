"""Equivalent fixed-graph SGConv using sparse matrix multiplication.

Only normalized, fixed adjacency is cached; trainable embeddings and their
gradients are recomputed on every call. Original model parameters are retained.
"""
import types
import torch
import torch.nn.functional as F
from torch_geometric.nn import SGConv
from torch_geometric.nn.conv.gcn_conv import gcn_norm

def forward(self, x, edge_index, edge_weight=None):
    assert not self.cached and self.flow == 'source_to_target'
    assert edge_weight is None or not edge_weight.requires_grad
    key=(x.shape[0], x.dtype, x.device, edge_index.data_ptr(),
         None if edge_weight is None else edge_weight.data_ptr())
    if getattr(self,'_wc_key',None) != key:
        e,w=gcn_norm(edge_index,edge_weight,x.shape[0],False,
                     self.add_self_loops,self.flow,dtype=x.dtype)
        self._wc_adj=torch.sparse_coo_tensor(e.flip(0),w,(x.shape[0],x.shape[0]),
                                             device=x.device).coalesce().to_sparse_csr()
        self._wc_key=key
    for _ in range(self.K):x=torch.sparse.mm(self._wc_adj,x)
    return self.lin(x)

def enable(model):
    count=0
    for layer in model.modules():
        if isinstance(layer,SGConv):
            layer.forward=types.MethodType(forward,layer);count+=1
    model.forward=types.MethodType(model_forward,model)
    return count

def model_forward(self,data):
    """Vectorize CopySlices and the final separable decoder, same equations."""
    x,pi=data.x,data.pert_idx
    n=len(pi);g=self.num_genes
    ids=torch.arange(g,device=x.device)
    # The published forward repeats gene-only embeddings for every cell before
    # the perturbation is added. Keep one copy, with exact multiplicities in BN.
    base=self.emb_trans(weighted_bn(self.bn_emb,self.gene_emb(ids),
                                    x.new_full((g,1),n)))
    pos=self.emb_pos(ids);other=pos
    for j,layer in enumerate(self.layers_emb_pos):
        pos=layer(pos,self.G_coexpress,self.G_coexpress_weight)
        # Original GEARS coexpression edges cover only the first graph in the
        # flattened batch. Later graphs have self loops only; preserve this.
        other=layer.lin(other)
        if j<len(self.layers_emb_pos)-1:pos=pos.relu();other=other.relu()
    proto=torch.cat([base+0.2*pos,base+0.2*other],0) if n>1 else base+0.2*pos
    weights=x.new_ones((len(proto),1))
    if n>1:weights[g:]=n-1
    for layer in self.emb_trans_v2.network:
        proto=weighted_bn(layer,proto,weights) if isinstance(layer,torch.nn.BatchNorm1d) else layer(proto)
    base=torch.cat([proto[:g][None],proto[g:][None].expand(n-1,-1,-1)],0) if n>1 else proto[None]
    pe=self.pert_emb(torch.arange(self.num_perts,device=x.device))
    for j,layer in enumerate(self.sim_layers):
        pe=layer(pe,self.G_sim,self.G_sim_weight)
        if j<self.num_layers-1:pe=pe.relu()
    active=[];sums=[]
    for j,indices in enumerate(pi):
        values=[pe[i] for i in indices if i!=-1]
        if values:
            active.append(j);total=values[0]
            for value in values[1:]:total=total+value
            sums.append(total)
    shifts=base.new_zeros((n,base.shape[-1]))
    if sums:
        inp=torch.stack(sums*2 if len(sums)==1 else sums)
        fused=self.pert_fuse(inp)[:len(sums)]
        shifts=shifts.index_add(0,torch.tensor(active,device=x.device),fused)
    base=self.bn_pert_base((base+shifts[:,None,:]).reshape(n*g,-1))
    rec=self.recovery_w(self.transform(base)).reshape(n,g,-1)
    single=(rec*self.indv_w1[:,:,0]).sum(2)+self.indv_b1[:,0]
    cross=self.cross_gene_state(single)
    # Sum over hidden dimensions algebraically, without repeating each cell's
    # cross-gene vector 5,006 times or constructing a concatenated dense tensor.
    result=single*self.indv_w2[0,:,0]+cross@self.indv_w2[0,:,1:].T+self.indv_b2
    return result+x.reshape(n,g)

def weighted_bn(layer,x,weight):
    if not layer.training:
        return F.batch_norm(x,layer.running_mean,layer.running_var,layer.weight,
                            layer.bias,False,0.,layer.eps)
    count=weight.sum();mean=(x*weight).sum(0)/count
    variance=((x-mean).square()*weight).sum(0)/count
    if layer.track_running_stats:
        with torch.no_grad():
            layer.num_batches_tracked.add_(1)
            momentum=layer.momentum if layer.momentum is not None else 1./float(layer.num_batches_tracked)
            layer.running_mean.lerp_(mean.detach(),momentum)
            layer.running_var.lerp_((variance*count/(count-1)).detach(),momentum)
    return (x-mean)*torch.rsqrt(variance+layer.eps)*layer.weight+layer.bias
