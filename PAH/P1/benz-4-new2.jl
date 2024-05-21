using FermiCG
using PyCall
using Plots
using LinearAlgebra
using Printf


pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("fcidump_4mer");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
ints = InCoreInts(ecore,h,g)

rdm1 = zeros(size(ints.h1))


na = 12
nb = 12

clusters_in    = [(1:6),(7:12),(13:18),(19:24)]
init_fspace = [(3,3),(3,3),(3,3),(3,3)]



# define clusters
clusters = [Cluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)

e_cmf, U, Da, Db  = FermiCG.cmf_oo(ints, clusters, init_fspace, rdm1,rdm1,
                                        max_iter_oo=20, verbose=0, gconv=1e-7, method="bfgs");
ints = FermiCG.orbital_rotation(ints,U)

max_roots = 100
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=0, max_roots=max_roots,
        init_fspace=init_fspace, rdm1a=Da, rdm1b=Db);

#cluster_bases = FermiCG.compute_cluster_est_basis(ints, clusters, Da, Db, thresh_schmidt=5e-5, init_fspace=init_fspace)

#
# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db, verbose=0);



nroots = 5

ref_fock = FermiCG.FockConfig(init_fspace)
ci_vector = FermiCG.ClusteredState(clusters, ref_fock, R=nroots)
ci_vector[ref_fock][ClusterConfig([2,1,1,1])] = [0,1,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,2,1,1])] = [0,0,1,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,2,1])] = [0,0,0,1,0]
ci_vector[ref_fock][ClusterConfig([1,1,1,2])] = [0,0,0,0,1]

thresh_list = [0.005,0.002,0.001,0.0007,0.0005,0.0003,0.0001]


for thresh_cipsi in thresh_list
    e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
    			    thresh_cipsi=thresh_cipsi, # Threshold for adding to P-space
    			    thresh_foi=1e-5,    # Threshold for keeping terms when defining FOIS    
    			    thresh_asci=0.00,     # Threshold of P-space configs to search from
    			    max_iter=10,
    			    matvec=3);

    e2a = FermiCG.compute_batched_pt2(v0, cluster_ops, clustered_ham, thresh_foi=1e-6)

    println()
    println("	*======TPSCI results======*")
    @printf("TCI Thresh: %8.6f  Dim:%8d\n",thresh_cipsi,size(v0)[1])
    println()
    @printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
    for r in 1:nroots
        @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2a[r] + ecore)
    end
    
    #global ci_vector = v0


    display(v0)
    end
end

