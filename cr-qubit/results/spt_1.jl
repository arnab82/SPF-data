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

v = FermiCG.BSTstate(clusters, FockConfig(init_fspace), cluster_bases)

xspace  = FermiCG.build_compressed_1st_order_state(v, cluster_ops, clustered_ham, nbody=3, thresh=1e-3)
xspace = FermiCG.compress(xspace, thresh=1e-3)
display(size(xspace))

FermiCG.nonorth_add!(v, xspace)
v = FermiCG.BSTstate(v,R=6)
FermiCG.randomize!(v)
FermiCG.orthonormalize!(v)

e_ci, v = FermiCG.ci_solve(v, cluster_ops, clustered_ham)

using Printf
display(v)
for ei in 1:length(e_ci)
    @printf(" %5i %12.8f %12.8f eV\n", ei, e_ci[ei]+ints.h0, (e_ci[ei]-e_ci[1])*27.21165)
end

 e_var, v_var = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,
                                        max_iter    = 20,
                                        nbody       = 4,
                                        H0          = "Hcmf",
                                        thresh_var  = 1e-2,
                                        thresh_foi  = 1e-4,
                                        thresh_pt   = 1e-3,
                                        ci_conv     = 1e-5,
                                        do_pt       = true,                          
                                        resolve_ss  = false,
                                        tol_tucker  = 1e-4)


@save "/home/arnab22/SPF-data/cr-qubit/data_spt1_1e2.jld2" e_var v_var
v1 = deepcopy(v_var);
display(v1,root=1)

 e2, v2 = FermiCG.block_sparse_tucker(v1, cluster_ops, clustered_ham,
                                        max_iter    = 20,
                                        nbody       = 4,
                                        H0          = "Hcmf",
                                        thresh_var  = 8e-3,
                                        thresh_foi  = 5e-5,
                                        thresh_pt   = 1e-3,
                                        ci_conv     = 1e-5,
                                        do_pt       = true,                          
                                        resolve_ss  = false,
                                        tol_tucker  = 1e-4)
@save "/home/arnab22/SPF-data/cr-qubit/data_spt1_5e3.jld2" e2 v2
display(v2, root=3, thresh=.1)
@time e_pt = FermiCG.compute_pt2_energy(v2, cluster_ops, clustered_ham, thresh_foi=1e-5,verbose=2,prescreen=true,compress_twice=true)

 e3, v3 = FermiCG.block_sparse_tucker(v2, cluster_ops, clustered_ham,
                                        max_iter    = 20,
                                        nbody       = 4,
                                        H0          = "Hcmf",
                                        thresh_var  = 6e-3,
                                        thresh_foi  = 2e-5,
                                        thresh_pt   = 8e-4,
                                        ci_conv     = 1e-5,
                                        do_pt       = true,                          
                                        resolve_ss  = false,
                                        tol_tucker  = 1e-4)

@save "/home/arnab22/SPF-data/cr-qubit/data_spt1_1e3.jld2" v3 e3
#@time e2 = FermiCG.compute_pt2_energy(v3, cluster_ops, clustered_ham, thresh_foi=1e-5)
@time e_pt = FermiCG.compute_pt2_energy(v3, cluster_ops, clustered_ham, thresh_foi=1e-5,verbose=2,prescreen=true,compress_twice=true)
