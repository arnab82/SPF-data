using QCBase
using Printf
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2

@load "cmf_data_2.jld2"
display(clusters)
display(init_fspace)
@load "cluster_bases_2.jld2"
@load "spt_var_2.jld2"
M=100
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3,5], FockConfig(init_fspace), max_roots=M, verbose=1);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

e_ci, v4 = FermiCG.ci_solve(v3, cluster_ops, clustered_ham)
#@time e2 = FermiCG.compute_pt2_energy2(v4, cluster_ops, clustered_ham, thresh_foi=1e-5,verbose=2,prescreen=true,compress_twice=true)
@time e2 = FermiCG.compute_pt2_energy(v4, cluster_ops, clustered_ham, thresh_foi=1e-5,verbose=2,prescreen=true,compress_twice=true)
