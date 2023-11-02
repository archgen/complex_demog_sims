#: Name: Matthew Williams
#: Date: 03:March:2023
# This snakemake script will run the WGS complex Eurasian demogrphic model from Williams et al. (2023) manuscript.

# Import python modules
import pandas as pd
import yaml

print("Config is: ", config)
print('\n' + '  ' * 80)

# --- Wildcards --- #
sim_model_name = config["simulation_model_name"]
iterations = [i for i in range(1,int(config["N_replicates"])+1)]
eig_files = ["snp", "ind", "geno"]
datatypes = ["sim_aDNA", "ts_branch"]
empir_aDNA_prefix_list = config["empir_aDNA_prefix"]


sample_sheet = pd.read_csv(config["sample_sheet"])
targetPops = list()
for i in range(sum(sample_sheet.qpAdmRotate == "target")):
    pop = ["{}".format(sample_sheet[sample_sheet.qpAdmRotate == "target"].Pop.values[i])]
    pop_str = str(pop).replace('[','').replace(']','').replace('\'','').replace('\"','')
    targetPops.append(pop_str.replace('[','').replace(']','').replace('\'','').replace('\"',''))

print('\n' + '-' * 80)
print("The analysis target populations are:\n")
for pop in targetPops:
    print(f"\t- {pop}")
print('\n' + '-' * 80)


rule all:
    input:
        # branch
        expand("out/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_admixtureF3_SummaryTable__{simName}.txt",  simName=sim_model_name),
        expand("out/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_FST_SummaryTable__{simName}.txt", simName=sim_model_name),
        expand("out/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_qpAdmReplacement_SummaryTable__{simName}.txt", simName=sim_model_name),
        expand("out/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_qpAdmRotation_SummaryTable__{simName}__{target}.txt", simName=sim_model_name, target=targetPops),
        # aDNA
        expand("out/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST_SummaryTable__{simName}.txt", simName=sim_model_name, empir_aDNA_prefix=empir_aDNA_prefix_list),
        expand("out/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3_SummaryTable__{simName}.txt", simName=sim_model_name, empir_aDNA_prefix=empir_aDNA_prefix_list),
        expand("out/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement_SummaryTable__{simName}.txt", simName=sim_model_name, empir_aDNA_prefix=empir_aDNA_prefix_list),
        expand("out/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation_SummaryTable__{simName}_{target}.txt", simName=sim_model_name, empir_aDNA_prefix=empir_aDNA_prefix_list, target=targetPops),

    shell:
        """
        echo -e "\033[48;5;45m\033[1;5;33m\033[3m\n      ☆\033[0m\033[31m\033[1mR\033[33mu\033[32mn\033[36m C\033[35mom\033[34mp\033[31ml\033[33met\033[0m\033[31m!\033[0m\033[48;5;45m\033[1;5;33m\033[3m ☆\n\033[0m"
        """


# This rule will generate a sim aDNA FST summary results from all the replicates.
rule ts_mu_eig_aDNA_FST_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST__{simName}_rep_{iteration}.rds", simName=sim_model_name, iteration=iterations, empir_aDNA_prefix=empir_aDNA_prefix_list)
    params:
        script = config["ts_mu_eig_aDNA_FST_Summary_Output"],
        data_dir = "out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST/"
    output:
        FST_sim_aDNA_Table = "out/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST_SummaryTable__{simName}.txt"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST_SummaryTable__{simName}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --in_dir {params.data_dir} \
        --out {output.FST_sim_aDNA_Table} "       


# --- Compute FST values from the simulated aDNA data --- #
# # This rule will generate an pairwise FST values for each individual replicate from the simulated aDNA data.
rule ts_mu_eig_aDNA_FST:
    input:
        sim_geno_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.geno",
        sim_snp_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.snp",
        sim_ind_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.ind"
    params:
        script = config["ts_mu_eig_aDNA_FST"],
        samples = config["sample_sheet"],
        sim_model = config["simulation_model_name"],
        datatype = "sim_aDNA"
    output:
        FST_aDNA = "out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST__{simName}_rep_{iteration}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST/ts_mu_eig_aDNA_{empir_aDNA_prefix}_FST__{simName}_rep_{iteration}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --sim_geno {input.sim_geno_aDNA} \
        --sample_sheet {params.samples} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --replicate {wildcards.iteration} \
        --out {output.FST_aDNA}"


# This rule will generate a summary table of the admixture F3 results from all the replicates and targets.
rule ts_mu_eig_aDNA_admixtureF3_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3__{simName}_rep_{iteration}__{target}.rds", simName=sim_model_name, iteration=iterations, target=targetPops, empir_aDNA_prefix=empir_aDNA_prefix_list)
    params:
        script = config["ts_mu_eig_aDNA_admixtureF3_SummaryTable"],
        workDIR = "out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3/"
    output:
        admixF3_Table = "out/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3_SummaryTable__{simName}.txt"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3_SummaryTable__{simName}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --in_dir {params.workDIR} \
        --out {output.admixF3_Table} "


# #--- Compute admixture F3 values from the simulated aDNA data --- #
# This rule will generate an admixture F3 analysis for each individual replicate from the simulated aDNA data
# on each of the target populations.
rule ts_mu_eig_aDNA_admixtureF3:
    input:
        sim_geno_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.geno",
        sim_snp_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.snp",
        sim_ind_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.ind",
    params:
        script = config["ts_mu_eig_aDNA_admixtureF3"],
        sample_sheet = config["sample_sheet"],
        datatype = "sim_aDNA",
        sim_model = config["simulation_model_name"]
    output:
        admixtureF3_aDNA = "out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3__{simName}_rep_{iteration}__{target}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3/ts_mu_eig_aDNA_{empir_aDNA_prefix}_admixtureF3__{simName}_rep_{iteration}__{target}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --sim_geno {input.sim_geno_aDNA} \
        --samples {params.sample_sheet} \
        --target {wildcards.target} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --replicate {wildcards.iteration} \
        --out {output.admixtureF3_aDNA}"


# # This rule will generate a qpAdm Replacement summary results table from all the replicates and all target populations.
rule ts_mu_eig_aDNA_qpAdmReplacement_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement_Table__{simName}_rep_{iteration}__{target}.rds", simName=sim_model_name, iteration=iterations, target=targetPops, empir_aDNA_prefix=empir_aDNA_prefix_list)
    params:
        script = config["ts_mu_eig_aDNA_qpAdmReplacement_SummaryTable"],
        sim_model = config["simulation_model_name"],
        datatype = "sim_aDNA",
        data_dir = "out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement/"
    output:
        qpAdmTable = "out/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement_SummaryTable__{simName}.txt"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement_SummaryTable__{simName}.benchmark.txt",
    shell:
        "Rscript {params.script} \
        --in_dir {params.data_dir} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --out {output.qpAdmTable}"


# --- Compute 2 Source Replacement qpAdm analyses from simulated aDNA data --- #
# # This rule will generate a qpAdm Replacement analysis for each individual replicate from the simulated aDNA data for all target populations.
rule ts_mu_eig_aDNA_qpAdmReplacement:
    input:
        sim_geno_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.geno",
        sim_snp_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.snp",
        sim_ind_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.ind"
    params:
        script = config["ts_mu_eig_aDNA_qpAdmReplacement"],
        sample_sheet = config["sample_sheet"]
    output:
        qpAdm2SourceFixed_table = "out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement_Table__{simName}_rep_{iteration}__{target}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmReplacement_Table__{simName}_rep_{iteration}__{target}.benchmark.txt",
    shell:
        "Rscript {params.script} \
        --sim_geno {input.sim_geno_aDNA} \
        --samples {params.sample_sheet} \
        --target {wildcards.target} \
        --iteration {wildcards.iteration} \
        --out {output.qpAdm2SourceFixed_table}"


#--- Process sim aDNA qpAdm Rotation Analyses --- #
# This rule will generate a qpAdm rotation summary results from all the replicates and targets.
rule ts_mu_eig_aDNA_qpAdmRotation_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation_{simName}_rep_{iteration}__{target}.rds", simName=sim_model_name, iteration=iterations, target=targetPops, empir_aDNA_prefix=empir_aDNA_prefix_list)
    params:
        data_dir = "out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation/",
        script_summary = config["ts_mu_eig_aDNA_qpAdmRotation_SummaryTable"],
        script_processing = config["qpAdmRotation_tables_processing"],
        sim_model = config["simulation_model_name"],
        datatype = "sim_aDNA",
        target=targetPops
    output:
        qpAdmTable = "out/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation_SummaryTable__{simName}_{target}.txt"
    benchmark:
        "benchmarks/{simName}/admixtools2/SummaryOutput/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation_SummaryTable__{simName}_{target}.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script_summary} \
        --in_dir {params.data_dir} \
        --target_pop {wildcards.target} \
        --processing_scr {params.script_processing} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --out {output.qpAdmTable} "


# #--- Compute qpAdm analyses from the simulated aDNA data --- #
# This rule will generate a qpAdm rotation analysis for each individual replicate from the simulated aDNA data 
# on the target population of interest.
rule ts_mu_eig_aDNA_qpAdmRotation:
    input:
        sim_geno_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.geno",
        sim_snp_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.snp",
        sim_ind_aDNA = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.ind",
    params:
        script = config["ts_mu_eig_aDNA_qpAdmRotate_loop"],
        samples = config["sample_sheet"]
    output:
        qpAdm_rotate_sim_aDNA = "out/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation_{simName}_rep_{iteration}__{target}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_mu_eig_aDNA/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation/ts_mu_eig_aDNA_{empir_aDNA_prefix}_qpAdmRotation_{simName}_rep_{iteration}__{target}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --sim_geno {input.sim_geno_aDNA} \
        --sim_ind {input.sim_ind_aDNA} \
        --sample_sheet {params.samples} \
        --target {wildcards.target} \
        --out {output.qpAdm_rotate_sim_aDNA}"


# --- Processing aDNA Conditions on mutated tree-sequence --- #
rule ts_mu_eig_aDNA_processing:
    input:
        sim_geno = "out/{simName}/ts_mu_eig/ts_mu_eig_{simName}_rep_{iteration}.geno",
        sim_snp = "out/{simName}/ts_mu_eig/ts_mu_eig_{simName}_rep_{iteration}.snp",
        sim_ind = "out/{simName}/ts_mu_eig/ts_mu_eig_{simName}_rep_{iteration}.ind",
    params:
        script = config["ts_mu_eig_aDNA_processing_lc"],
        script_func = config["aDNA_conditions_functions_lc"],
        samples = config["sample_sheet"],
        aDNA_data_dir = config["empir_aDNA_dir"],
        empir_prefix = "{empir_aDNA_prefix}",
        out_prefix = "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}" 
    output:
        "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.geno",
        "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.snp",
        "out/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.ind"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/ts_mu_eig_aDNA_{empir_aDNA_prefix}/ts_mu_eig_aDNA_{empir_aDNA_prefix}_{simName}_rep_{iteration}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --sim_geno {input.sim_geno} \
        --sim_snp {input.sim_snp} \
        --sim_ind {input.sim_ind} \
        --empir_data_prefix {params.empir_prefix} \
        --aDNA_data_dir {params.aDNA_data_dir} \
        --sample_sheet {params.samples} \
        --script_func {params.script_func} \
        --out {params.out_prefix} "


# --- Output Eigenstrat Formatted Data from the mutated tree-sequence --- #
rule ts_mu_eigenstrat:
    input:
        ts = "out/{simName}/ts_mu/ts_mu_{simName}_rep_{iteration}.ts",
    params:
        script = config["msp_ts_mu_eig_output"],
        samples = config["sample_sheet"],
        sim_model = config["simulation_model_name"],
        outDIR = "out/{simName}/ts_mu_eig/"
    output:
        "out/{simName}/ts_mu_eig/ts_mu_eig_{simName}_rep_{iteration}.geno",
        "out/{simName}/ts_mu_eig/ts_mu_eig_{simName}_rep_{iteration}.snp",
        "out/{simName}/ts_mu_eig/ts_mu_eig_{simName}_rep_{iteration}.ind"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/ts_mu_eig/ts_mu_eig_{simName}_rep_{iteration}.benchmark.txt"
    shell:
        "python3 {params.script} \
        --simName {params.sim_model} \
        --tree_sequence {input.ts} \
        --samples {params.samples} \
        --iteration {wildcards.iteration} \
        --out {params.outDIR}"

# --- Simulate mutations on msprime tree-sequence --- #
rule msp_ts_mutate:
    input:
        ts = "out/{simName}/ts/ts_{simName}_rep_{iteration}.ts"
    params:
        script = config["msp_ts_mutate"],
        mu = config["mutation_rate"],
        DTWF_time = config["DTWF_generations"]
    output:
        ts_mu = "out/{simName}/ts_mu/ts_mu_{simName}_rep_{iteration}.ts"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/ts_mu/ts_mu_{simName}_rep_{iteration}.benchmark.txt"
    shell: 
        "python3 {params.script} \
        --tree_sequence {input.ts} \
        --mutation_rate {params.mu} \
        --WFgens {params.DTWF_time} \
        --iteration {wildcards.iteration} \
        --out {output.ts_mu} "


# This rule will generate an admixture F3 summary results from all the replicates and target popualtions.
rule ts_branch_admixtureF3_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_branch/ts_branch_admixtureF3/ts_branch_admixtureF3__{simName}_rep_{iteration}__{target}.rds", simName=sim_model_name, iteration=iterations, target=targetPops)
    params:
        script = config["ts_branch_admixtureF3_SummaryTable"],
        sim_model = config["simulation_model_name"],
        replicates = config["N_replicates"],
        data_dir = "out/{simName}/admixtools2/ts_branch/ts_branch_admixtureF3/"
    output:
        admixF3_Table = "out/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_admixtureF3_SummaryTable__{simName}.txt"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_admixtureF3_SummaryTable__{simName}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --in_dir {params.data_dir} \
        --out {output.admixF3_Table} "


# #--- Compute admixture F3 from F2 matrix branch lengths --- #
# This rule will generate admixture F3 analyses for each individual replicate from the F2 matrix generated on the simulated branch lengths
# for each target population of interest.
rule ts_branch_admixtureF3:
    input:
        f2_matrix = "out/{simName}/admixtools2/ts_branch/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.rds"
    params:
        script = config["ts_branch_admixtureF3"],
        samples = config["sample_sheet"],
        datatype = "ts_branch",
        sim_model = config["simulation_model_name"]
    output:
        admixtureF3_F2_branch = "out/{simName}/admixtools2/ts_branch/ts_branch_admixtureF3/ts_branch_admixtureF3__{simName}_rep_{iteration}__{target}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch/ts_branch_admixtureF3/ts_branch_admixtureF3__{simName}_rep_{iteration}__{target}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --f2_matrix {input.f2_matrix} \
        --samples {params.samples} \
        --target {wildcards.target} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --replicate {wildcards.iteration} \
        --out {output.admixtureF3_F2_branch}"


# This rule will generate a qpAdm Replacement summary results table from all the replicates and all target populations.
rule ts_branch_qpAdmReplacement_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_branch/ts_branch_qpAdmReplacement_Table/ts_branch_qpAdmReplacement_Table__{simName}_rep_{iteration}__{target}.rds", simName=sim_model_name, iteration=iterations, target=targetPops)
    params:
        script = config["ts_branch_qpAdmReplacement_SummaryTable"],
        sim_model = config["simulation_model_name"],
        datatype = "ts_branch_F2",
        data_dir = "out/{simName}/admixtools2/ts_branch/ts_branch_qpAdmReplacement_Table/"
    output:
        qpAdmTable = "out/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_qpAdmReplacement_SummaryTable__{simName}.txt"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_qpAdmReplacement_SummaryTable__{simName}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --in_dir {params.data_dir} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --out {output.qpAdmTable}"


# #--- Compute qpAdm Replacement analyses from F2 matrix branch lengths --- #
# This rule will generate a qpAdm Replacement analysis for each individual replicate from the F2 matrix generated on the simulated branch lengths.
rule ts_branch_qpAdmReplacement:
    input:
        script = config["ts_branch_qpAdmReplacement"],
        f2_matrix = "out/{simName}/admixtools2/ts_branch/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.rds"
    params:
        samples = config["sample_sheet"]
    output:
        qpAdmReplacement_table = "out/{simName}/admixtools2/ts_branch/ts_branch_qpAdmReplacement_Table/ts_branch_qpAdmReplacement_Table__{simName}_rep_{iteration}__{target}.rds"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch/ts_branch_qpAdmReplacement_Table/ts_branch_qpAdmReplacement_Table__{simName}_rep_{iteration}__{target}.benchmark.txt",
    shell:
        "Rscript {input.script} \
        --f2_matrix {input.f2_matrix} \
        --samples {params.samples} \
        --target {wildcards.target} \
        --iteration {wildcards.iteration} \
        --out {output.qpAdmReplacement_table}"


#--- Process Branch qpAdm Rotation Analyses --- #
# This rule will generate a qpAdm rotation summary results from all the replicates and targets.
rule ts_branch_qpAdmRotation_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_branch/ts_branch_qpAdmRotation/ts_branch_qpAdmRotation__{simName}_rep_{iteration}__{target}.rds", simName=sim_model_name, iteration=iterations, target=targetPops)
    params:
        data_dir = "out/{simName}/admixtools2/ts_branch/ts_branch_qpAdmRotation/",
        script_summary = config["ts_branch_qpAdmRotation_SummaryTable"],
        script_processing = config["qpAdmRotation_tables_processing"],
        sim_model = config["simulation_model_name"],
        datatype = "ts_branch_F2",
        target=targetPops
    output:
        qpAdmTable = "out/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_qpAdmRotation_SummaryTable__{simName}__{target}.txt"
    benchmark:
        "benchmarks/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_qpAdmRotation_SummaryTable__{simName}__{target}.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script_summary} \
        --in_dir {params.data_dir} \
        --target_pop {wildcards.target} \
        --processing_scr {params.script_processing} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --out {output.qpAdmTable} "


#--- Compute qpAdm analyses from F2 matrix branch lengths --- #
# This rule will generate a qpAdm rotation analysis for each individual replicate from the F2 matrix generated on the simulated branch lengths
# on the target popualtion of interest.
rule ts_branch_qpAdmRotation:
    input:
        f2_matrix = "out/{simName}/admixtools2/ts_branch/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.rds"
    params:
        script = config["ts_branch_qpAdmRotation"],
        samples = config["sample_sheet"], 
    output:
        qpAdm_rotate_F2_branch = "out/{simName}/admixtools2/ts_branch/ts_branch_qpAdmRotation/ts_branch_qpAdmRotation__{simName}_rep_{iteration}__{target}.rds"
    wildcard_constraints:
       iteration="\d+"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch/ts_branch_qpAdmRotation/ts_branch_qpAdmRotation__{simName}_rep_{iteration}__{target}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script} \
        --f2_matrix {input.f2_matrix} \
        --samples {params.samples} \
        --target {wildcards.target} \
        --out {output.qpAdm_rotate_F2_branch}"


# This rule will generate an F2 matrix for each individual replicate ts from the simulated branch lengths.
rule ts_branch_F2_matrix:
    input:
        f2_table = "tmp/{simName}/admixtools2/ts_branch/ts_branch_F2_DT/ts_branch_F2_DT__{simName}_rep_{iteration}.txt"
    params:
        script = config["ts_branch_F2_matrix"]
    output:
        f2_matrix = "out/{simName}/admixtools2/ts_branch/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.rds"
    wildcard_constraints:
       iteration="\d+"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch/ts_branch_F2_matrix/ts_branch_F2_matrix__{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "Rscript {params.script} \
        --f2_table {input.f2_table} \
        --out {output.f2_matrix}"


#--- Compute tables of F2 statistics from branch lengths --- #
# This rule will generate an F2 table for each individual replicate ts from the simulated branch lengths.
rule ts_branch_F2_DT:
    input:
        expand("out/{simName}/ts/ts_{simName}_rep_{iteration}.ts", simName=sim_model_name, iteration=iterations)
    params:
        script = config["ts_branch_F2_DT"],
        window_size = config["window_length_Fstats"],
        ts = "out/{simName}/ts/ts_{simName}_rep_{iteration}.ts"
    output:
        f2_table = "tmp/{simName}/admixtools2/ts_branch/ts_branch_F2_DT/ts_branch_F2_DT__{simName}_rep_{iteration}.txt"
    wildcard_constraints:
       iteration="\d+"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch/ts_branch_F2_DT/ts_branch_F2_DT__{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "python3 {params.script} \
        --tree_sequence {params.ts} \
        --window_size {params.window_size} \
        --out {output.f2_table}"


# This rule will generate FST summary results from all the replicates.
rule ts_branch_FST_Summary_Output:
    input:
        expand("out/{simName}/admixtools2/ts_branch/ts_branch_FST_DT/ts_branch_FST_DT__{simName}_rep_{iteration}.rds", simName=sim_model_name, iteration=iterations)
    params:
        script = config["ts_branch_FST_Summary_Output"],
        sim_model = config["simulation_model_name"],
        datatype = "ts_branch_F2",
        data_dir = "out/{simName}/admixtools2/ts_branch/ts_branch_FST_DT/"
    output:
        FST_sim_aDNA_Table = "out/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_FST_SummaryTable__{simName}.txt"
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/admixtools2/SummaryOutput/ts_branch/ts_branch_SummaryOutput_FST_SummaryTable__{simName}.benchmark.txt"
    shell:
        "Rscript {params.script} \
        --in_dir {params.data_dir} \
        --data_type {params.datatype} \
        --model_name {params.sim_model} \
        --out {output.FST_sim_aDNA_Table}"


#--- Compute matrix of FST statistics from branch lengths --- #
# This rule will generate an FST matrix for each individual replicate ts from the simulated branch lengths.
rule ts_branch_FST_DT:
    input:
        expand("out/{simName}/ts/ts_{simName}_rep_{iteration}.ts", simName=sim_model_name, iteration=iterations)
    params:
        script = config["ts_branch_FST_DT"],
        simName = config["simulation_model_name"],
        ts = "out/{simName}/ts/ts_{simName}_rep_{iteration}.ts",
    output:
        fst_table = "out/{simName}/admixtools2/ts_branch/ts_branch_FST_DT/ts_branch_FST_DT__{simName}_rep_{iteration}.rds"
    wildcard_constraints:
       iteration="\d+"
    benchmark:
        "benchmarks/{simName}/admixtools2/ts_branch/ts_branch_FST_DT/ts_branch_FST_DT__{simName}_rep_{iteration}.benchmark.txt"
    conda:
        "msprime-env.yaml"
    shell:
        "python3 {params.script} \
        --tree_sequence {params.ts} \
        --simName {params.simName} \
        --iteration {wildcards.iteration} \
        --out {output.fst_table}"


# --- Simulate msprime tree-sequence --- #
rule msp_ts_simulate:
    input:
        script = config["msp_ts_simulate"],
        demes = config["demes_yaml"],
        samples = config["sample_sheet"],
    params:
        simName = config["simulation_model_name"],
        DTWF_time = config["DTWF_generations"]
    output:
        ts = "out/{simName}/ts/ts_{simName}_rep_{iteration}.ts" 
    wildcard_constraints:
       iteration="\d+"
    conda:
        "msprime-env.yaml"
    benchmark:
        "benchmarks/{simName}/ts/ts_{simName}_rep_{iteration}.benchmark.txt"
    shell:
        "python3 {input.script} \
        --simName {params.simName} \
        --demes {input.demes} \
        --samples {input.samples} \
        --WFgens {params.DTWF_time} \
        --iteration {wildcards.iteration} \
        --out {output.ts} "


# ---- ENV DETAILS --- #
rule conda_r_python_env_libpaths_info:
    conda:
        "msprime-env.yaml"
    shell:
        r"""
        echo "Conda environment: $CONDA_PREFIX"
        Rscript -e 'cat("R library paths:", .libPaths(), sep = "\n")'
        python -c 'import sys; print("Python library paths:", *sys.path, sep = "\n")'
        """
