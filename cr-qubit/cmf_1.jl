using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers
#load integrals from disk
h0 = npzread("clustered1_integrals_h0.npy")
h1 = npzread("clustered1_integrals_h1.npy")
h2 = npzread("clustered1_integrals_h2.npy")
ints = InCoreInts(h0, h1, h2)
ints_original = deepcopy(ints);

println(" Integrals have the following sizes: h0= ", size(h0), " h1= ", size(h1), " h2= ", size(h2))


# Define clusters - this should probably be done in the python notebook, and then just read in here
clusters_in = [
	       (1:6), # Benzene 1
	       (7:12), # Benzene 1
	       (13:18), # Benzene 1
	       (19:24), # Benzene 1
	       (25:34),   # metal
	       ]

clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]

init_fspace = [
	       (3,3),
	       (3,3),
	       (3,3),
	       (3,3),
	       (5,5)];
display(clusters)
display(init_fspace)


rdm1 = zeros(size(ints.h1))

e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(rdm1, rdm1), verbose=0, gconv=1e-6, method="bfgs", sequential=true);

ints = RDM.orbital_rotation(ints_original,U);
d1=RDM.orbital_rotation(d1,U);
@save "cmf_data_1.jld2" ints clusters init_fspace d1