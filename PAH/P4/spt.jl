using FermiCG
using LinearAlgebra
using Printf
using JLD2
@load "data_cmf_qdphenyl.jld2"

max_roots = 100
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1_n, [3,3,3,3], FermiCG.FockConfig(init_fspace), max_roots=max_roots, verbose=1);


# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
nroots = 17

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1_n.a, d1_n.b)
ci_vector = FermiCG.BSTstate(clusters, FockConfig(init_fspace), cluster_bases,R=nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.TuckerConfig((1:1,1:1,1:1,1:1))] =
    FermiCG.Tucker(tuple([zeros(Float64, 1, 1,1,1) for _ in 1:nroots]...))
FermiCG.add_single_excitons!(ci_vector,FermiCG.FockConfig(init_fspace),nroots)
FermiCG.add_spin_flip_states!(ci_vector,FermiCG.FockConfig(init_fspace),1)
FermiCG.add_1electron_transfers!(ci_vector,FermiCG.FockConfig(init_fspace),1)
FermiCG.eye!(ci_vector)

e_ci, v = FermiCG.ci_solve(ci_vector, cluster_ops, clustered_ham,conv_thresh = 6e-5);
thresh_list = [0.01,0.008,0.006]

thresh_spf=0.01
#for thresh_spf in thresh_list
   e_var, v_var = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 8e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               verbose     = 1)
    @save "/home/arnab22/SPF-data/P4/data_spt.jld2" e_var v_var cluster_bases

thresh_spf=0.009
e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 8e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               verbose     = 1)
    @save "/home/arnab22/SPF-data/P4/data_spt_2.jld2" e_var v_var cluster_bases
thresh_spf=0.0085

e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/100,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 8e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               verbose     = 1)
    @save "/home/arnab22/SPF-data/P4/data_spt_3.jld2" e_var v_var cluster_bases
thresh_spf=0.008

e_var, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 20,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/100,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 8e-5,
                                               do_pt       = false,
                                               resolve_ss  = false,
                                               verbose     = 1)
    @save "/home/arnab22/SPF-data/P4/data_spt_3.jld2" e_var v_var cluster_bases
