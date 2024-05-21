using FermiCG, NPZ
using JLD2
@load "cmf_data_2.jld2"
display(clusters)
display(init_fspace)
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3,1,4], FockConfig(init_fspace), max_roots=20, verbose=1);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);
@save "cluster_bases_2.jld2" cluster_bases
v = FermiCG.BSTstate(clusters, FockConfig(init_fspace), cluster_bases, R=3)

FermiCG.add_single_excitons!(v, FockConfig(init_fspace))
FermiCG.randomize!(v)
FermiCG.orthonormalize!(v)
e_ci, v = FermiCG.ci_solve(v, cluster_ops, clustered_ham);
