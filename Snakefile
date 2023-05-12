#should be run from Intros-CH-AlphaDelta

# Cluster-specific settings for resources required by any rule. This file
# provides default resources for all rules and allows users to specify resources
# per rule by name. An important resource for the Hutch cluster is the requested
# "partition". Jobs submitted to the "restart" partition will start running
# almost immediately, but they may also be killed at any moment when someone
# else needs those resources. This is analogous to the spot resources on AWS.
#cluster-config: "cluster.json"

configfile: "config.yaml"

rule all_pie_delta:
    input:
        ["results/Delta/trees/tree_data{runn}.csv".format(runn=runn) for runn in range(0,10)] #runs 0-9

rule all_pie_alpha:
    input:
        ["results/Alpha/trees/tree_data{runn}.csv".format(runn=runn) for runn in range(0,10)] #runs 0-9

#rule two_pie_delta:
#    input:
#        ["results/Delta/trees/tree_data{runn}.csv".format(runn=runn) for runn in range(8,10)] #runs 0-9


rule all_clus_alpha:
    input:
        ["results/Alpha/clusters/liberalClusters-{runn}.csv".format(runn=runn) for runn in range(0,10)] #runs 0-9

rule all_clus_delta:
    input:
        ["results/Delta/clusters/liberalClusters-{runn}.csv".format(runn=runn) for runn in range(0,10)] #runs 0-9


#folder format name is too different so need a helper function
def _alpha_or_delta_folder(wildcards):
    if wildcards.var_type == "Delta":
        return f"21A.Delta-swiss_{wildcards.run}-2021-07-31"
    else:
        return f"20I.Alpha.V1-swiss_{wildcards.run}-2021-03-31"

rule run_pies:
    message: "Run alpha or delta builds in pie plot"
    input:
        metadata = "../../ncov_2021_AlphDelt/data/metadata.tsv", #from the builds
    params:
        run_folder = _alpha_or_delta_folder
    output:
        "results/{var_type}/trees/tree_data{run}.csv"
    log:
        "logs/{var_type}/pie_{run}.txt"
    shell:
        """
        python3 scripts/tree_pie_plot.py \
            --run-folder {params.run_folder} \
            --variant-type {wildcards.var_type} \
            --nextstrain-metadata {input.metadata} 2>&1 | tee {log}
        """

rule run_clus_intros:
    message: "Run alpha or delta builds to determine lib & cons clusters"
    input:
        "results/{var_type}/out_data/pie_slices{run}.json"
    params:
        results_folder = "results/"
    output:
        lib = "results/{var_type}/clusters/liberalClusters-{run}.csv",
        con = "results/{var_type}/clusters/conservClusters-{run}.csv",
    log:
        "logs/{var_type}/clus_{run}.txt" 
    shell:
        """
        python3 scripts/analyze_slices.py \
            --run-folder {params.results_folder} \
            --variant-type {wildcards.var_type} \
            --run-number {wildcards.run} 2>&1 | tee {log}
        """

#rule run_pie_alphas:
#    message: "Run alpha builds in pie plot"
#    input:
#        metadata = "../../ncov_2021_AlphDelt/data/metadata.tsv", #from the builds
#        run_folder = "20I.Alpha.V1-swiss_{run}-2021-03-31"
#    output:
#        "results/Alpha/trees/tree_data{run}.csv"
#    log:
#        "logs/{var_type}/pie_{run}.txt"
#    shell:
#        """
#        python3 scripts/tree_pie_plot.py \
#            --run-folder {input.run_folder} \
#            --variant-type {var_type} \
#            --nextstrain-metadata {input.metadata} 2>&1 | tee {log}
#        """

#rule add_labels:
#    message: "Remove extraneous colorings for main build and move frequencies"
#    input:
#        #auspice_json = rules.incorporate_travel_history.output.auspice_json,
#        auspice_json = rules.include_hcov19_prefix.output.auspice_json,
#        tree = rules.refine.output.tree,
#        clades = rules.clades.output.clade_data,
#        mutations = rules.ancestral.output.node_data
#    output:
#        auspice_json = "results/{build_name}/ncov_with_accessions_and_travel_branches_and_labels.json",
#    log:
#        "logs/add_labels_{build_name}.txt"
#    conda: config["conda_environment"]
#    shell:
#        """
#        python3 scripts/add_labels.py \
#            --input {input.auspice_json} \
#            --tree {input.tree} \
#            --mutations {input.mutations} \
#            --clades {input.clades} \
#            --output {output.auspice_json} 2>&1 | tee {log}
#        """