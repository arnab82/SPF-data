using FermiCG
using LinearAlgebra
using Printf
using JLD2


na = 12
nb = 12

clusters_in    = [(1:6),(7:12),(13:18),(19:24)]
init_fspace = [(3,3),(3,3),(3,3),(3,3)]



# define clusters
clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)


@load "data_cmf_qdphenyl.jld2"
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
   
e_ci, v_ci = FermiCG.ci_solve(v, cluster_ops, clustered_ham)
display(e_ci)
println("*************************************************************************************************************************************************************************************")

v = FermiCG.BSTstate(v,R=1)
xspace  = FermiCG.build_compressed_1st_order_state(v, cluster_ops, clustered_ham, nbody=4, thresh=1e-3)
xspace = FermiCG.compress(xspace, thresh=1e-3)

FermiCG.nonorth_add!(v, xspace)
v = FermiCG.BSTstate(v,R=6)
FermiCG.randomize!(v)
FermiCG.orthonormalize!(v)
e_ci, v = FermiCG.ci_solve(v, cluster_ops, clustered_ham)
println("*************************************************************************************************************************************************************************************")

e_pt, v_pt = FermiCG.do_fois_pt2(v, cluster_ops, clustered_ham, thresh_foi=1e-3, max_iter=50, tol=1e-8)
println("*************************************************************************************************************************************************************************************")
    

e_var, v_var = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,
                                               max_iter    = 50,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = 1e-2,
                                               thresh_foi  = 1e-3,
                                               thresh_pt   = 1e-3,
                                               ci_conv     = 1e-5,
                                               do_pt       = true,
                                               resolve_ss  = false,
                                               tol_tucker  = 1e-4,
                                               solver      = "davidson")

