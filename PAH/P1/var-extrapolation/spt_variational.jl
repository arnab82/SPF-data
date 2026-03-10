using FermiCG
using LinearAlgebra
using Printf
using JLD2
@load "data_benz_p1_cmf.jld2"
max_roots = 60
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1_n, [3,3,3,3], FermiCG.FockConfig(init_fspace), max_roots=max_roots, verbose=1);

# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
nroots = 1

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1_n.a, d1_n.b)

v = FermiCG.BSTstate(clusters, FockConfig(init_fspace), cluster_bases)
thresh_spf=0.01
e_var1, v_var = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var1 v_var cluster_bases

thresh_spf=0.008
e_var2, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var2 v_var cluster_bases



thresh_spf=0.006
e_var3, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var3 v_var cluster_bases




thresh_spf=0.004
e_var4, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var4 v_var cluster_bases

thresh_spf=0.002
e_var4, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var4 v_var cluster_bases


thresh_spf=0.001
e_var4, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var4 v_var cluster_bases
