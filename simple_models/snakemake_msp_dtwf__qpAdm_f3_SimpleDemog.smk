#: Name: Matthew Williams
#: Date: 20:Feb:2023
# This snakemake script will run the simple demographic Models from the Williams et al. (2023) manuscript. 

# Import python packages
import pandas as pd
import yaml

# Print the input config file to screen
print("Config is: ", config)

# --- Wildcards --- #
sim_model_name = config["simulation_model_name"]
iterations = [i for i in range(1,int(config["N_replicates"])+1)]

# Rule all - output files at the end of the pipeline dependency structure.
rule all:
    input:
        expand("out/{simName}/ts/ts_{simName}_rep_{iteration}.ts", simName=sim_model_name, iteration=iterations),
        expand("out/{simName}/admixtools2/ts_branch_SummaryOutput_qpAdmRotation/ts_branch_qpAdmrotation_SummaryTable__{simName}.txt", simName=sim_model_name),
        expand("out/{simName}/admixtools2/ts_branch_SummaryOutput_qpAdmRotation_2source/ts_branch_qpAdmrotation_2source_SummaryTable__{simName}.txt", simName=sim_model_name),
        expand("out/{simName}/admixtools2/ts_branch_SummaryOutput_FST/ts_branch_FST_SummaryTable__{simName}.rds", simName=sim_model_name, iteration=iterations),
        expand("out/{simName}/admixtools2/ts_branch_SummaryOutput_admixtureF3/ts_branch_admixtureF3_SummaryTable__{simName}.rds", simName=sim_model_name)


# This rule will generate an admixture F3 summary results from all the replicates.
rule ts_branch_admixtureF3_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_branch_admixtureF3/ts_branch_admixtureF3__{simName}_rep_{iteration}.rds", simName=sim_model_name, iteration=iterations)
    params:
        script = config["ts_branch_admixF3_SummaryTable"],
        sim_model = config["simulation_model_name"],
        replicates = config["N_replicates"],
        datatype = "ts_branch",
        workDIR = "out/{simName}/admixtools2/ts_branch_admixtureF3/"
    output:
        admixF3_Table = "out/{simName}/admixtools2/ts_branch_SummaryOutput_admixtureF3/ts_branch_admixtureF3_SummaryTable__{simName}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_SummaryOutput_admixtureF3/ts_branch_admixtureF3_SummaryTable__{simName}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --in_dir {params.workDIR} \
        --data_type {params.datatype} \
        --replicates {params.replicates} \
        --model_name {params.sim_model} \
        --out {output.admixF3_Table} "


#--- Compute qpAdm analyses from F2 matrix branch lengths --- #
# This rule will generate a qpAdm rotation analysis for each individual replicate from the F2 matrix generated on the simulated branch lengths.
rule ts_branch_admixtureF3:
    input:
        f2_matrix = "out/{simName}/admixtools2/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.rds",
    params:
        script = config["ts_branch_admixF3"]
    output:
        admixF3_branch = "out/{simName}/admixtools2/ts_branch_admixtureF3/ts_branch_admixtureF3__{simName}_rep_{iteration}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_admixtureF3/ts_branch_admixtureF3__{simName}_rep_{iteration}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --f2_matrix {input.f2_matrix} \
        --replicate {wildcards.iteration} \
        --out {output.admixF3_branch}"

#--- Process Branch qpAdm Analyses & generate Figures / Tables --- #
# This rule will generate a qpAdm rotation summary results from all the replicates.
rule ts_branch_qpAdmRotation_2source_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_branch_qpAdmrotation_2source/ts_branch_qpAdmrotation_2source__{simName}_rep_{iteration}.rds", simName=sim_model_name, iteration=iterations)
    params:
        fileDIR = "out/{simName}/admixtools2/ts_branch_qpAdmrotation_2source/",
        script = config["ts_branch_qpAdmRotation_SummaryTable_2source"],
        RscriptDIR = config["R_script_DIR"],
        sim_model = config["simulation_model_name"],
        datatype = "ts_branch_F2"
    output:
        qpAdmTable = "out/{simName}/admixtools2/ts_branch_SummaryOutput_qpAdmRotation_2source/ts_branch_qpAdmrotation_2source_SummaryTable__{simName}.txt"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_SummaryOutput_qpAdmRotation_2source/ts_branch_qpAdmrotation_2source_SummaryTable__{simName}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script} \
        --in_dir {params.fileDIR} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --RscriptDIR {params.RscriptDIR} \
        --out {output.qpAdmTable} "


#--- Compute qpAdm analyses from F2 matrix branch lengths --- #
# This rule will generate a qpAdm rotation analysis for each individual replicate from the F2 matrix generated on the simulated branch lengths.
rule ts_branch_qpAdmRotation_2source:
    input:
        f2_matrix = "out/{simName}/admixtools2/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.rds"
    params:
        script = config["ts_branch_qpAdmRotation_2source"]
    output:
        qpAdm_rotate_F2_branch = "out/{simName}/admixtools2/ts_branch_qpAdmrotation_2source/ts_branch_qpAdmrotation_2source__{simName}_rep_{iteration}.rds"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_qpAdmrotation/ts_branch_qpAdmrotation_2source__{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script} \
        --f2_matrix {input.f2_matrix} \
        --out {output.qpAdm_rotate_F2_branch}"


#--- Process Branch qpAdm Analyses & generate Figures / Tables --- #
# This rule will generate a qpAdm rotation summary results from all the replicates.
rule ts_branch_qpAdmRotation_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_branch_qpAdmrotation/ts_branch_qpAdmrotation__{simName}_rep_{iteration}.rds", simName=sim_model_name, iteration=iterations)
    params:
        fileDIR = "out/{simName}/admixtools2/ts_branch_qpAdmrotation/",
        script = config["ts_branch_qpAdmRotation_SummaryTable"],
        RscriptDIR = config["R_script_DIR"],
        sim_model = config["simulation_model_name"],
        datatype = "ts_branch_F2"
    output:
        qpAdmTable = "out/{simName}/admixtools2/ts_branch_SummaryOutput_qpAdmRotation/ts_branch_qpAdmrotation_SummaryTable__{simName}.txt"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_SummaryOutput_qpAdmRotation/ts_branch_qpAdmrotation_SummaryTable__{simName}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script} \
        --in_dir {params.fileDIR} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --RscriptDIR {params.RscriptDIR} \
        --out {output.qpAdmTable} "


#--- Compute qpAdm analyses from F2 matrix branch lengths --- #
# This rule will generate a qpAdm rotation analysis for each individual replicate from the F2 matrix generated on the simulated branch lengths.
rule ts_branch_qpAdmRotation:
    input:
        f2_matrix = "out/{simName}/admixtools2/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.rds"
    params:
        script = config["ts_branch_qpAdmRotation"]
    output:
        qpAdm_rotate_F2_branch = "out/{simName}/admixtools2/ts_branch_qpAdmrotation/ts_branch_qpAdmrotation__{simName}_rep_{iteration}.rds"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_qpAdmrotation/ts_branch_qpAdmrotation__{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script} \
        --f2_matrix {input.f2_matrix} \
        --out {output.qpAdm_rotate_F2_branch}"


#--- Compute matrix of F2 statistics from branch lengths --- #
# This rule will generate an F2 matrix for each individual replicate ts from the simulated branch lengths.
rule ts_branch_F2_matrix:
    input:
        f2_table = "out/{simName}/admixtools2/ts_branch_F2_table/ts_branch_F2_table__{simName}_rep_{iteration}.txt"
    params:
        script = config["ts_branch_F2_matrix"]
    output:
        f2_matrix = "out/{simName}/admixtools2/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.rds"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script} \
        --f2_table {input.f2_table} \
        --out {output.f2_matrix}"


#--- Compute tables of F2 statistics from branch lengths --- #
# This rule will generate an F2 table for each individual replicate ts from the simulated branch lengths.
rule ts_branch_F2:
    input:
        expand("out/{simName}/ts/ts_{simName}_rep_{iteration}.ts", simName=sim_model_name, iteration=iterations)
    params:
        script = config["ts_branch_F2"],
        window_size = config["window_length_Fstats"],
        ts = "out/{simName}/ts/ts_{simName}_rep_{iteration}.ts"
    output:
        f2_table = "out/{simName}/admixtools2/ts_branch_F2_table/ts_branch_F2_table__{simName}_rep_{iteration}.txt"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_F2_table/ts_branch_F2_table__{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "python3 {params.script} \
        --tree_sequence {params.ts} \
        --window_size {params.window_size} \
        --out {output.f2_table}"


# This rule will generate a branch average FST summary results from all the replicates.
rule ts_branch_FST_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_branch_FST/ts_branch_FST_table__{simName}_rep_{iteration}.rds", simName=sim_model_name, iteration=iterations)
    params:
        script = config["ts_branch_FST_Summary_Output"],
        sim_model = config["simulation_model_name"],
        inDIR = "out/{simName}/admixtools2/ts_branch_FST/"
    output:
        FST_sim_aDNA_Table = "out/{simName}/admixtools2/ts_branch_SummaryOutput_FST/ts_branch_FST_SummaryTable__{simName}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_SummaryOutput_FST/ts_branch_FST_SummaryTable__{simName}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --in_dir {params.inDIR} \
        --out {output.FST_sim_aDNA_Table}"


#--- Compute table of FST statistics from branch lengths --- #
# This rule will generate a pairwise FST table for each individual replicate ts from the simulated branch lengths.
rule ts_branch_FST:
    input:
        expand("out/{simName}/ts/ts_{simName}_rep_{iteration}.ts", simName=sim_model_name, iteration=iterations)
    params:
        script = config["ts_branch_FST"],
        simName = config["simulation_model_name"],
        ts = "out/{simName}/ts/ts_{simName}_rep_{iteration}.ts",
    output:
        fst_table = "out/{simName}/admixtools2/ts_branch_FST/ts_branch_FST_table__{simName}_rep_{iteration}.rds"
    wildcard_constraints:
       iteration="\d+"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch_FST/ts_branch_FST_table__{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "python3 {params.script} \
        --tree_sequence {params.ts} \
        --simName {params.simName} \
        --iteration {wildcards.iteration} \
        --out {output.fst_table}"

# --- Simulate msprime tree-sequence --- #
rule ts_simulate:
    params:
        script = config["msp_ts_simulate"],
        simName = config["simulation_model_name"],
        seq_length = config["chromosome_length"],
        recom_rate = config["recombination_rate"],
        DTWFgens = config["DTWF_generations"],
        MRCA_gen = config["MRCA_split_generation"],
        N_diploid_samples = config["N_diploid_samples"],
        Ne = config["diploid_Ne"],
        genTime = config["generation_time"]
    output:
        ts = "out/{simName}/ts/ts_{simName}_rep_{iteration}.ts"
    wildcard_constraints:
       iteration="\d+"
    benchmark:
        "benchmarks/{simName}/ts/ts_{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "python3 {params.script} \
        --simName {params.simName} \
        --length {params.seq_length} \
        --recom {params.recom_rate} \
        --iteration {wildcards.iteration} \
        --DTWFgens {params.DTWFgens} \
        --MRCAgen {params.MRCA_gen} \
        --genTime {params.genTime} \
        --Ne {params.Ne} \
        --Nsamples {params.N_diploid_samples} \
        --out {output.ts}"
