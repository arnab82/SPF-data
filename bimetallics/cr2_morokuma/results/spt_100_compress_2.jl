using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "data_cmf_29_cr2.jld2"
@load "data_bst_3e3.jld2"
M = 100 

init_fspace = FockConfig([(6, 3), (4,4), (3,6)])

#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/cluster_ops_100.jld2" cluster_ops
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4
display(vbst)
v_bst_3e4=FermiCG.compress(vbst, thresh=3e-4)
@time e2_bst = FermiCG.compute_pt2_energy(v_bst_3e4, cluster_ops, clustered_ham, thresh_foi=1e-6,prescreen=true,compress_twice=true,verbose=2);
v_bst_25e4=FermiCG.compress(vbst, thresh=2.5e-4)
@time e2_bst = FermiCG.compute_pt2_energy(v_bst_25e4, cluster_ops, clustered_ham, thresh_foi=1e-6,prescreen=true,compress_twice=true,verbose=2);

@time e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-6,prescreen=true,compress_twice=true,verbose=2);
