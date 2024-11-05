using FermiCG
using LinearAlgebra
using Printf
using JLD2
@load "data_benz-6_cmf.jld2"

max_roots = 100
# Build Cluster basis
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1_n, [3,3,3,3,3,3], FermiCG.FockConfig(init_fspace), max_roots=max_roots, verbose=1);
@load "data_spt.jld2"

# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
nroots = 25

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1_n.a, d1_n.b)

#e_ci, v_var = FermiCG.ci_solve(ci_vector, cluster_ops, clustered_ham,conv_thresh = 5e-4);
thresh_list = [0.01,0.009]

thresh_spf=0.01

    
    println("**********************************************************************")
e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 1,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/3,
                                               ci_conv     = 1e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               tol_tucker  = 1e-4,
                                       	       verbose     = 1)
@save "/home/arnab22/SPF-data/P5/data_spt.jld2" e_var v_var cluster_bases
e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 1,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/3,
                                               ci_conv     = 1e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               tol_tucker  = 1e-4,
                                               verbose     = 1)
@save "/home/arnab22/SPF-data/P5/data_spt.jld2" e_var v_var cluster_bases
