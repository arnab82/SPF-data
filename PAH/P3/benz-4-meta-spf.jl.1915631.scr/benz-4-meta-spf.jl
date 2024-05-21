using FermiCG
using PyCall
using LinearAlgebra
using Printf
using ClusterMeanField  
using ActiveSpaceSolvers
using RDM
pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("benz-4-meta-fcidump");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
ints = InCoreInts(ecore,h,g)

#rdm1 = zeros(size(ints.h1))
rdm1=RDM1(n_orb(ints))

na = 12
nb = 12

clusters_in    = [(1:6),(7:12),(13:18),(19:24)]
init_fspace = [(3,3),(3,3),(3,3),(3,3)]



# define clusters
clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)


ansatze=[FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3)]
e_cmf, U, d1_n = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace, ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-8, 
                           tol_d1=1e-10,
                           tol_ci=1e-11,
                           verbose=4, 
                           zero_intra_rots = true,
                           sequential=true)
ints = FermiCG.orbital_rotation(ints,U)

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

e_var, v_var = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,
                                               max_iter    = 20,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = 1e-2,
                                               thresh_foi  = 1e-3,
                                               thresh_pt   = sqrt(1e-5),
                                               ci_conv     = 1e-5,
                                               do_pt       = true,
                                               resolve_ss  = true,
                                               tol_tucker  = 1e-4, 
                                               verbose     = 1)
