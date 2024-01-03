using FermiCG
using PyCall
using LinearAlgebra
using Printf
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using JLD2

pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("benz-6-fcidump");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
ints = InCoreInts(ecore,h,g)

#rdm1 = zeros(size(ints.h1))
rdm1=RDM1(n_orb(ints))
na = 18
nb = 18

clusters_in    = [(1:6),(7:12),(13:18),(19:24),(25:30),(31:36)]
init_fspace = [(3,3),(3,3),(3,3),(3,3),(3,3),(3,3)]



# define clusters
clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)

ansatze=[FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3),FCIAnsatz(6,3,3)]
e_cmf, U, d1_n = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace, ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-8,
                           tol_d1=1e-10,
                           tol_ci=1e-11,
                           verbose=4,
                           zero_intra_rots = true,
                           sequential=true)
ints = FermiCG.orbital_rotation(ints,U)
@save "data_benz-6_cmf.jld2" clusters init_fspace ints d1_n e_cmf ansatze U ecore
