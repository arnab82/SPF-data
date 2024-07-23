using FermiCG
using LinearAlgebra
using Printf
using JLD2
@load "data_benz-p3-cmf.jld2"
@load "data_spt.jld2"
max_roots = 100
# Build Cluster basis
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1_n, [3,3,3,3], FermiCG.FockConfig(init_fspace), max_roots=max_roots, verbose=1);
thresh_list = [0.008,0.007,0.006]

   # Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

   # Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
nroots = 17

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1_n.a, d1_n.b)
display(v_var)
v6=FermiCG.compress(v_var, thresh=6e-3)
display(v6)
@time e2 = FermiCG.compute_pt2_energy(v6, cluster_ops, clustered_ham, thresh_foi=1e-4,verbose=2,prescreen=true,compress_twice=true);

v5=FermiCG.compress(v_var, thresh=5e-3)
display(v5)
@time e2 = FermiCG.compute_pt2_energy(v5, cluster_ops, clustered_ham, thresh_foi=1e-4,verbose=2,prescreen=true,compress_twice=true);

v4=FermiCG.compress(v_var, thresh=4e-3)
display(v4)
@time e2 = FermiCG.compute_pt2_energy(v4, cluster_ops, clustered_ham, thresh_foi=1e-4,verbose=2,prescreen=true,compress_twice=true);

v3=FermiCG.compress(v_var, thresh=3e-3)
display(v3)
@time e2 = FermiCG.compute_pt2_energy(v3, cluster_ops, clustered_ham, thresh_foi=1e-4,verbose=2,prescreen=true,compress_twice=true);

@time e2 = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=1e-4,verbose=2,prescreen=true,compress_twice=true);
