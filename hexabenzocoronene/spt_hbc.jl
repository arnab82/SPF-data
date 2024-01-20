using FermiCG
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals   
using JLD2 
using Printf
using NPZ 
@load "/Users/ayush/workspace/cMF_data/pi_conjugated_systems/data_cmf_pi_cg_4.jld2"

display(clusters)
display(init_fspace)
max_roots = 40

#
# form Cluster data
cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=0,
                                                   max_roots=max_roots,
                                                   init_fspace=init_fspace,
                                                   rdm1a=d1.a, rdm1b=d1.b)

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

 
v = FermiCG.BSTstate(clusters, FockConfig(init_fspace), cluster_bases)


xspace  = FermiCG.build_compressed_1st_order_state(v, cluster_ops, clustered_ham, nbody=4, thresh=1e-3)
xspace = FermiCG.compress(xspace, thresh=1e-3)
display(size(xspace))

FermiCG.nonorth_add!(v, xspace)
v = FermiCG.BSTstate(v,R=3)
FermiCG.randomize!(v)
FermiCG.orthonormalize!(v)

e_ci, v = FermiCG.ci_solve(v, cluster_ops, clustered_ham)

display(v)
for ei in 1:length(e_ci)
    @printf(" %5i %12.8f %12.8f eV\n", ei, e_ci[ei]+ints.h0, (e_ci[ei]-e_ci[1])*27.21165)
end


e_var, v_var = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,
                                            max_iter    = 20,
                                            nbody       = 4,
                                            H0          = "Hcmf",
                                            thresh_var  = 1e-1,
                                            thresh_foi  = 1e-3,
                                            thresh_pt   = 1e-3,
                                            ci_conv     = 1e-5,
                                            do_pt       = true,                          
                                            resolve_ss  = false,
                                            tol_tucker  = 1e-4)

display(v_var)
@time e2 = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=1e-5)