using QCBase
using Printf
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load  "cmf_diis.jld2"
@load "spt_50_var.jld2"
M = 75

display(clusters)
display(init_fspace)
ints=deepcopy(ints_cmf)
ref_fspace = FockConfig(init_fspace)
ecore = ints.h0
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
#cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
display(v_var)
@load "cluster_ops_50.jld2"
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

#display(v_var)
@time ept2_thresh_2 = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=2e-4,verbose=2,prescreen   = true,compress_twice = true)
@time ept2_thresh_3 = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=3e-4,verbose=2,prescreen   = true,compress_twice = true)

@time ept2_thresh_4 = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=4e-4,verbose=2,prescreen   = true,compress_twice = true)
v_spt_2=FermiCG.compress(v_var, thresh=2e-3)
display(v_spt_2)
@time ept2_2 = FermiCG.compute_pt2_energy(v_spt_2, cluster_ops, clustered_ham, thresh_foi=1e-4,verbose=2,prescreen   = true,compress_twice = true)


v_spt_3=FermiCG.compress(v_var, thresh=3e-3)
display(v_spt_3)
@time ept2_3 = FermiCG.compute_pt2_energy(v_spt_3, cluster_ops, clustered_ham, thresh_foi=1e-4,verbose=2,prescreen   = true,compress_twice = true)

v_spt_4=FermiCG.compress(v_var, thresh=4e-3)
display(v_spt_4)
@time ept2_4 = FermiCG.compute_pt2_energy(v_spt_4, cluster_ops, clustered_ham, thresh_foi=1e-4,verbose=2,prescreen   = true,compress_twice = true)
