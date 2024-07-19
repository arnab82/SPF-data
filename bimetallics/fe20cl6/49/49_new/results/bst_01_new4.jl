using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

@load "data_cmf_49_new.jld2"

M = 100


init_fspace =  [(5, 0), (4, 4), (0, 5), (12, 12), (0, 0)]
init_fspace = FermiCG.FockConfig(init_fspace)
@load "data_spt_new.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6



bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham, thresh_var=8e-3, thresh_foi=0.8e-4,thresh_pt=2.5e-3,thresh_spin=1.5e-3,ci_conv=4e-5,ci_max_iter=500,do_pt=false);
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/49_new/data_spt_new.jld2" clusters init_fspace ints cluster_bases v

#FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-6)

bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham, thresh_var=8e-3, thresh_foi=0.6e-4,thresh_pt=1e-3,thresh_spin=1e-3,ci_conv=4e-5,ci_max_iter=500,do_pt=false);
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/49_new/data_spt_new.jld2" clusters init_fspace ints cluster_bases v

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-6)

