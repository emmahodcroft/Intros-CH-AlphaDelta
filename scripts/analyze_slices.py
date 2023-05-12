import json
import re
import os, sys
from collections import defaultdict
import pandas as pd
from Bio import Phylo
from datetime import datetime
from statistics import mean
from collections import Counter

# Create a dictionary to find nodes by name
def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

#calculate mean date from a list of datetime objects
def get_mean_date(dates):
    timestampDates = [x.timestamp() for x in dates]
    mn = sum(timestampDates) / len(dates)
    return datetime.fromtimestamp(mn)

#recursively print swiss parent nodes and parents - also number them
def printNodeAndChildren(tabs, node, clusterNumber):
    global parentsToShow
    #if child is not in parent list, we are done - just assign number
    if node not in parentsToShow:
        names[node].cladeNum = clusterNumber
        return
    print(tabs*"\t",node,parentList[node], clusterNumber)
    for no in parentList[node]: #try assigning to children
        names[no].cladeNum = clusterNumber 
    names[node].cladeNum = clusterNumber
    parentsToShow.remove(node)
    for c in parentList[node]:
        printNodeAndChildren(tabs+1, c, clusterNumber)

#######################


if __name__=="__main__":
    import argparse

    parser = parser = argparse.ArgumentParser(description='Analayze collapsed "pies" to find potential introductions',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run-folder', help="path to 'results' folder of pie collapsing run")
    parser.add_argument('--variant-type', help="Specify whether Alpha or Delta build")
    parser.add_argument('--run-number', help="Number of run to process")
    parser.add_argument('--print-nodes', action='store_true', help="set to print out visual representation of swiss nodes")
    args = parser.parse_args()

    build = args.variant_type
    runNum = args.run_number
    
    #create paths for input files
    pie_slices_file = f"{args.run_folder}/{build}/out_data/pie_slices{runNum}.json"
    meta_file = f"{args.run_folder}/{build}/trees/tree_data{runNum}.csv"
    tree_file = f"{args.run_folder}/{build}/trees/treeSubtree{runNum}.nwk"

    #check these exist before we proceed further!
    for fi in [pie_slices_file, meta_file, tree_file]:
        if not os.path.exists(fi):
            print("FILE NOT FOUND!")
            print(fi)
            sys.exit()

    ############
    #create output files
    #results already exists so create other folders
    clusters_folder = f"{args.run_folder}/{build}/clusters"
    if not os.path.isdir(clusters_folder):
        print(f"\nCouldn't find folder {clusters_folder} - creating it.\n")
        os.mkdir(clusters_folder)
    liberal_file = f"{clusters_folder}/liberalClusters-{runNum}.csv"
    conserv_file = f"{clusters_folder}/conservClusters-{runNum}.csv"
    #for those with all swiss dates listed
    liberal_file_dates = f"{clusters_folder}/liberalClusters_dates-{runNum}.csv"
    conserv_file_dates = f"{clusters_folder}/conservClusters_dates-{runNum}.csv"
    #log file
    log_file = f"{clusters_folder}/log-{runNum}.txt"

    log_text = ""

    ############
    #read in data

    #read in the pie slice data
    with open(pie_slices_file) as f:
        data = json.load(f)

    #read in the tree, add onto it the metdata:
    meta = pd.read_csv(meta_file, sep=',', index_col=0)
    tree = Phylo.read(tree_file, "newick") 

    #be able to find nodes in tree by name
    names = lookup_by_names(tree)

    #link children to parents & vice versa
    for node in tree.find_clades(order='preorder'):
        if not node.is_terminal():
            node.internal_children = []
        #if child in node:
            node.direct_children = []
        for child in node.clades:
            child.parent = node
            if child.is_terminal():
                node.direct_children.append(child)
            else:
                node.internal_children.append(child)

    #attach data to tree
    for node in tree.find_clades(order='postorder'):
        if node.name in meta.index:
            node.country = meta.loc[node.name, "country"]
            node.division =meta.loc[node.name, "division"]
            node.date = meta.loc[node.name, "date"]
        else:
            node.country = ''
            node.division = ''
            node.date = ''

    ####
    # Start processing

    swissNodes = len(data["Switzerland"])
    #Nodes that contain swiss sequences once collapsed
    log_text += f"Nodes that contain Swiss sequences once collapsed: {swissNodes}\n"
    # Alpha Apr 22 - 1038
    # Delta Apr 22 - 1347


    #But some of these have parents that also contain Swiss Seqs:
    numSwissPar = 0
    for k,v in data["Switzerland"].items():
        if 'parent' in v:
            numSwissPar+=1
    log_text += f"\nNumber of parents that also have Swiss seqs: {numSwissPar}\n"
    # Alpha Apr 22 - 655
    # Delta Apr 22 - 892



    swissOnlyNum = len(data["Switzerland"])-numSwissPar
    log_text += f"\nNumber Swiss-only with no Swiss seq in parents: {swissOnlyNum}\n"
    # Alpha Apr 22 - 383
    # Delta Apr 22 - 455 are Swiss-only with no Swiss seq in parents

    # Lets see how big they are
    #currently don't look at this in command-line version as comparison hard.
    sizeSwissOnly = {}
    sizes = []
    for k,v in data["Switzerland"].items():
        if 'parent' not in v:
            sizeSwissOnly[k] = len(v['sequences'])
            sizes.append(len(v['sequences']))

    sizeWithParents = {}
    sizesP = []
    for k,v in data["Switzerland"].items():
        if 'parent' in v:
            sizeWithParents[k] = len(v['sequences'])
            sizesP.append(len(v['sequences']))

    # swiss only have no children in parents, because if they did they
    # would not be swiss only

    # what about those where the parent is the same?
    parentList = defaultdict(list)
    p = re.compile("NODE_[0-9]+")
    for k,v in data["Switzerland"].items():
        if 'parent' in v:
            #grep out parent name
            m = p.findall(v["parent"])
            parent = m[0]
            parentList[parent].append(k)

    log_text += f"\nThere are {len(parentList)} parents for the {numSwissPar} clusters that have Swiss parents\n"
    #Alpha Apr 22 - there are 345 parents for the 655 clusters that have swiss parents
    #Delta Apr 22 - there are 379 parents for the 892 clusters that have swiss parents

    #lets find those which share parents
    sharedParents = {}
    for k,v in parentList.items():
        if len(v) > 1:
            sharedParents[k] = v

    # create child to parent dict to find parents easily
    childList = {}
    for k,v in parentList.items():
        for c in v:
            childList[c] = k


    # list of parents, to remove those we have already printed
    # and number the clades & children - CHECK THIS WORKED
    parentsToShow = list(parentList.keys())
    rootParents =[]
    clusterNumber = 0
    for k,v in parentList.items():
        node = k
        # if node has a parent, back up the tree
        while node in childList and childList[node] in parentsToShow:
            node = childList[node]
        if node in parentsToShow:
            rootParents.append(node)
            if(args.print_nodes):
                printNodeAndChildren(0, node, clusterNumber)
                print("")
            clusterNumber += 1

    lenRootPts = len(rootParents)
    lenNoSwissPts = len([x for x in rootParents if x in sizeSwissOnly.keys()])
    # check all of the rootParents are also nodes without Swiss parents (as expected)
    log_text += f"\nThere are {lenRootPts} rootParents and {lenNoSwissPts} nodes without Swiss parents (should be same)\n"


    ###COUNTING THE CLUSTERS!!

    #we want to print info in two ways - one counting conservatively (rootParents only once)
    # and one counting liberaly (each with swiss parents counted alone)

    #let's do liberal first, easier?
    liberalClusterInfo = defaultdict(dict)
    for k,v in data["Switzerland"].items():
        #record clade number
        cl = names[k].cladeNum if hasattr(names[k], "cladeNum") else pd.NA
        #find out if 'top' of clade
        headClusterNode = True if k in rootParents else pd.NA

        directChildren = v['sequences']
        numberDirectChildren = len(directChildren)

        #get information about Swiss children
        datesDirectChildren = [datetime.strptime(names[ch].date, "%Y-%m-%d") for ch in directChildren]
        meanDate = get_mean_date(datesDirectChildren).strftime("%Y-%m-%d")
        minDate = min(datesDirectChildren).strftime("%Y-%m-%d")
        maxDate = max(datesDirectChildren).strftime("%Y-%m-%d")
        allSwissDates = [d.strftime("%Y-%m-%d") for d in datesDirectChildren]
        #get canton stats
        canton_dict = dict(Counter([names[ch].division for ch in directChildren]))
        child_cantons = ', '.join(['{}: {}'.format(*i) for i in canton_dict.items()])

        #get info about other children
        all_children = names[k].direct_children 
        numberAllChildren = len(all_children)

        nonSwissChildren = [c for c in all_children if c.country != "Switzerland"]
        if nonSwissChildren:
            countriesNonSwissChildren_dict = dict(Counter([c.country for c in nonSwissChildren]))
            countriesNonSwissChildren = ', '.join(['{}: {}'.format(*i) for i in countriesNonSwissChildren_dict.items()])

            datesNonSwissChilds = [datetime.strptime(ch.date, "%Y-%m-%d") for ch in nonSwissChildren]
            meanDateNS = get_mean_date(datesNonSwissChilds).strftime("%Y-%m-%d")
            minDateNS = min(datesNonSwissChilds).strftime("%Y-%m-%d")
            maxDateNS = max(datesNonSwissChilds).strftime("%Y-%m-%d")

            #figure out earliest non-swiss Seq
            [c.country for c in nonSwissChildren if datetime.strptime(c.date, "%Y-%m-%d") == min(datesNonSwissChilds)]

        else:
            countriesNonSwissChildren = ""
            meanDateNS = ""
            minDateNS = ""
            maxDateNS = ""

        liberalClusterInfo[k] = {
            "Clade": cl,
            "HeadOfClade": headClusterNode,
            "NumberSwissChildren": numberDirectChildren,
            "TotalNumberChildren": numberAllChildren,

            "MeanDateSwissChildren": meanDate,
            "MinDateSwissChildren": minDate,
            "MaxDateSwissChildren": maxDate,
            "CantonsSwissChildren": child_cantons,

            "NonSwissChildrenCountries": countriesNonSwissChildren,
            "MeanDateNonSwissChildren": meanDateNS,
            "MinDateNonSwissChildren": minDateNS,
            "MaxDateNonSwissChildren": maxDateNS,

            "TotalDatesSwissChildren": allSwissDates,
        }

    #convert to dataframe
    liberalClusterDF = pd.DataFrame.from_dict(data=liberalClusterInfo, orient='index')
    liberalClusterDF.sort_values(by=["Clade","HeadOfClade"]).to_csv(liberal_file_dates, header=True)
    liberalClusDF_nodate = liberalClusterDF.drop("TotalDatesSwissChildren", axis=1)
    liberalClusDF_nodate.sort_values(by=["Clade","HeadOfClade"]).to_csv(liberal_file, header=True)

    #shortcut to print straight out (unordered)
    #(pd.DataFrame.from_dict(data=liberalClusterInfo, orient='index')
    #    .to_csv('liberalClusters.csv', header=True))

    conservativeClusterInfo = defaultdict(dict)
    #use the entries from before with no (null/NA) cluster number
    swissOnlyInfo = liberalClusterDF[liberalClusterDF.isnull().Clade]
    #ane entries from before who are heads of cluster
    clusterHeadInfo = liberalClusterDF[liberalClusterDF.notnull().HeadOfClade]

    conservClusterDF = pd.concat([swissOnlyInfo, clusterHeadInfo])
    conservClusterDF.sort_values(by=["Clade","HeadOfClade"]).to_csv(conserv_file_dates, header=True)
    conservClusterDF.drop("TotalDatesSwissChildren", axis=1, inplace=True)
    conservClusterDF.sort_values(by=["Clade","HeadOfClade"]).to_csv(conserv_file, header=True)

    with open(log_file, 'w') as f:
        f.write(log_text)

    print(f"Run {runNum} completed.")


