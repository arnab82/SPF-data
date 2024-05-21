using FermiCG
using LinearAlgebra
using Printf
using JLD2
@load "data_benz_p1_cmf.jld2"

max_roots = 100
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=0, max_roots=max_roots,
        init_fspace=init_fspace, rdm1a=d1_n.a, rdm1b=d1_n.b);


# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1_n.a, d1_n.b)
v = FermiCG.BSTstate(clusters, FockConfig(init_fspace), cluster_bases)
thresh_list = [0.005,0.002,0.001,0.0007,0.0005,0.0003,0.0001]


for thresh_spf in thresh_list
   e_var, v_var = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,
                                               max_iter    = 20,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = 1e-5,
                                               thresh_pt   = sqrt(1e-5),
                                               ci_conv     = 1e-5,
                                               do_pt       = true,
                                               resolve_ss  = true,
                                               tol_tucker  = 1e-6, 
                                               verbose     = 1)
    e_pt = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=1e-5)
    println("**********************************************************************")
    println("   *======SPF results======*")
    @printf("SPF Thresh: %8.6f  Dim:%8d\n",thresh_spf,size(v_var)[1])
    println()
    @printf("SPF %5s %12s %12s\n", "Root", "E(0)","E(2)")
    @printf("TCI %5s %12.8f %12.8f\n",1, e_var[1] + ecore,ecore+e_pt[1])
    println("**********************************************************************")
end
