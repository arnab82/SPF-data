using QCBase
using Printf
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2

@load "cmf_data_1.jld2"
display(clusters)
display(init_fspace)
@load "cluster_bases_1.jld2"
M=100
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3,5], FockConfig(init_fspace), max_roots=M, verbose=1);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

@load "spt_1_var.jld2"

 e3, v3 = FermiCG.block_sparse_tucker(v3, cluster_ops, clustered_ham,
                                        max_iter    = 20,
                                        nbody       = 4,
                                        H0          = "Hcmf",
                                        thresh_var  = 1e-2,
                                        thresh_foi  = 1e-5,
                                        thresh_pt   = 8e-4,
					thresh_spin = 8e-3,
					ci_max_iter =  100,
                                        ci_conv     = 1e-5,
                                        do_pt       = false,                          
                                        resolve_ss  = false,
                                        tol_tucker  = 1e-4)

@save "/home/arnab22/SPF-data/cr-qubit/data_spt1_final.jld2" v3 e3
#@time e2 = FermiCG.compute_pt2_energy(v3, cluster_ops, clustered_ham, thresh_foi=1e-5)
@time e_pt = FermiCG.compute_pt2_energy(v3, cluster_ops, clustered_ham, thresh_foi=5e-5,verbose=2,prescreen=true,compress_twice=true)
