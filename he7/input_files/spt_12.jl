using QCBase
using Printf
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


function run_spt_scan()
    for geom in reverse(12:12)
        f_out = @sprintf("%03i.out",geom)
        f_err = @sprintf("%03i.err",geom)

        #redirect_stdio(stdout=f_out, stderr=f_err) do
            file1 = "data_cmf_$(geom).jld2"
            println(" CMF File: ", file1)
            @eval @load $file1
            println(" Enuc: ", cmf_ints.h0)

            M = 100

            init_fspace = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
            ref_fspace = FockConfig(init_fspace)
            ecore = cmf_ints.h0

            cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(cmf_ints, clusters, d1, [1,1,1,1,1,1,1], ref_fspace, max_roots=M, verbose=4);

            clustered_ham = FermiCG.extract_ClusteredTerms(cmf_ints, clusters)
            cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, cmf_ints);

            FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, cmf_ints, d1.a, d1.b);
            @load "data_spt_012.jld2"
            @time e2 = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=1e-6,prescreen   = true,compress_twice = true)
    
            e_var .+= cmf_ints.h0
            eguess .+= cmf_ints.h0
            e2guess .+= cmf_ints.h0
            e2 .+= cmf_ints.h0

            print("*EGuess: ")
            [@printf("%12.8f ", i) for i in eguess]
            println()

            print("*E2Guess: ")
            [@printf("%12.8f ", i) for i in e2guess]
            println()

            print("*E0: ")
            [@printf("%12.8f ", i) for i in e_var]
            println()

            print("*E2: ")
            [@printf("%12.8f ", i) for i in e2]
            println()

            @save @sprintf("data_spt_%03i.jld2",geom) vguess v_var  eguess e_var e2 e2guess 
        #end
    end
end

run_spt_scan()


