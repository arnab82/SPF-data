using FermiCG
using LinearAlgebra
using Printf
using JLD2
@load "data_biphenyl_2mer_1d_cmf.jld2"

max_roots = 100
# Build Cluster basis
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1_n, [3,3,3,3], FermiCG.FockConfig(init_fspace), max_roots=max_roots, verbose=1);
@load "data_spt.jld2"

# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
nroots = 17

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1_n.a, d1_n.b)

thresh_list = [0.01,0.009]

thresh_spt=0.0085
thresh_spf=0.009
   e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 5,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spt,
                                               thresh_foi  = thresh_spf/45,
                                               thresh_pt   = thresh_spt/2.5,
                                               ci_conv     = 1e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
					       tol_tucker  = 5e-5,
                                               verbose     = 1)
@save "/home/arnab22/SPF-data/P2/data_spt.jld2" e_var v_var cluster_bases
   e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 5,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spt,
                                               thresh_foi  = thresh_spf/45,
                                               thresh_pt   = thresh_spt/2.5,
                                               ci_conv     = 1e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               tol_tucker  = 5e-5,
                                               verbose     = 1)
@save "/home/arnab22/SPF-data/P2/data_spt.jld2" e_var v_var cluster_bases
   e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 5,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spt,
                                               thresh_foi  = thresh_spf/45,
                                               thresh_pt   = thresh_spt/2.5,
                                               ci_conv     = 1e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               tol_tucker  = 5e-5,
                                               verbose     = 1)
@save "/home/arnab22/SPF-data/P2/data_spt.jld2" e_var v_var cluster_bases
   e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 5,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spt,
                                               thresh_foi  = thresh_spf/45,
                                               thresh_pt   = thresh_spt/2.5,
                                               ci_conv     = 1e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               tol_tucker  = 5e-5,
                                               verbose     = 1)
@save "/home/arnab22/SPF-data/P2/data_spt.jld2" e_var v_var cluster_bases
   e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 5,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spt,
                                               thresh_foi  = thresh_spf/45,
                                               thresh_pt   = thresh_spt/2.5,
                                               ci_conv     = 1e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               tol_tucker  = 5e-5,
                                               verbose     = 1)
@save "/home/arnab22/SPF-data/P2/data_spt.jld2" e_var v_var cluster_bases

@time e2 = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=1e-5, verbose=2,prescreen=true,compress_twice=true)

v1=FermiCG.compress(v_var,thresh=4e-3)
@time e2 = FermiCG.compute_pt2_energy(v1, cluster_ops, clustered_ham, thresh_foi=1e-5, verbose=2,prescreen=true,compress_twice=true)

v2=FermiCG.compress(v_var,thresh=5e-3)
@time e2 = FermiCG.compute_pt2_energy(v2, cluster_ops, clustered_ham, thresh_foi=1e-5, verbose=2,prescreen=true,compress_twice=true)

v3=FermiCG.compress(v_var,thresh=6e-3)
@time e2 = FermiCG.compute_pt2_energy(v3, cluster_ops, clustered_ham, thresh_foi=1e-5, verbose=2,prescreen=true,compress_twice=true)

v4=FermiCG.compress(v_var,thresh=7e-3)
@time e2 = FermiCG.compute_pt2_energy(v4, cluster_ops, clustered_ham, thresh_foi=1e-5, verbose=2,prescreen=true,compress_twice=true)
