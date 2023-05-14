import json
import re
from collections import defaultdict
import pandas as pd
from Bio import Phylo
from datetime import datetime
from statistics import mean
from collections import Counter

#This is same as the `_orig` version of the script, to be run locally and to
#work with the original file structure of early runs for the manuscript
#modified to produce output files that also contain all dates of Swiss sequences
#to help with plotting early cluster date distributions

#run in 2021_Delta

build = "Delta"

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

#read in the pie slice data
with open(f'out_data/{build}/pie_slices.json') as f:
  data = json.load(f)

#read in the tree, add onto it the metdata:
meta = pd.read_csv(f"treeFiles/{build}/tree_data.csv", sep=',', index_col=0)
tree = Phylo.read(f"treeFiles/{build}/treeSubtree.nwk", "newick") 

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

# take a look
data["Switzerland"]

#Nodes that contain swiss sequences once collapsed
len(data["Switzerland"])
# Alpha Apr 22 - 1038
# Delta Apr 22 - 1347

# Alpha Feb - 953
# Alpha Nov - 409
# Delta Nov -  120
# Europe july - 101

#But some of these have parents that also contain Swiss Seqs:

numSwissPar = 0
for k,v in data["Switzerland"].items():
    if 'parent' in v:
        numSwissPar+=1

numSwissPar
# Alpha Apr 22 - 655
# Delta Apr 22 - 892

# Alpha Feb - 608
# Alpha Nov - 333
# Delta Nov - 82
# Europe July -  39 have parents which also contain Swiss seqs

len(data["Switzerland"])-numSwissPar
# Alpha Apr 22 - 383
# Delta Apr 22 - 455

# Alpha Feb - 345
# Alpha Nov - 76
# Delta Nov - 38
# Europe July - 62 are Swiss-only with no Swiss seq in parents

# Lets see how big they are
sizeSwissOnly = {}
sizes = []
for k,v in data["Switzerland"].items():
    if 'parent' not in v:
        sizeSwissOnly[k] = len(v['sequences'])
        sizes.append(len(v['sequences']))
#Alpha Apr 22
#[1,1,1,1,1,4,4,1,1,2,24,17,6,31,1,2,9,5,2,5,4,1,1,6,39,1,1,2,1,160,1,55,1,3,1,1,3,2,27,1,1,2,1,2,27,1,1,13,1,6,3,1,2,1,1,17,1,1,1,3,5,1,1,3,1,61,1,2,1,4,5,1,3,1,1,19,1,33,1,1,6,1,2,5,1,1,1,1,1,1,1,2,1,1,1,1,9,1,2,1,1,1,3,1,3,4,1,1,15,1,2,1,1,2,2,2,30,1,3,1,1,1,2,3,1,1,8,2,1,1,24,7,2,2,4,4,1,10,1,14,3,2,1,1,1,5,30,1,1,1,1,1,2,5,1,20,3,13,6,3,2,1,28,1,2,1,3,1,1,1,1,2,1,6,3,1,1,23,10,10,1,2,1,13,4,7,1,9,1,3,1,11,1,1,1,11,2,1,1,1,1,2,1,1,3,1,3,3,1,1,3,1,4,32,2,1,2,8,1,1,1,2,1,68,2,9,2,2,8,1,5,1,2,8,6,1,1,1,1,2,5,1,1,1,3,64,1,1,8,4,1,1,1,1,4,1,1,2,1,11,4,1,1,4,2,2,1,20,2,1,6,1,40,3,1,3,1,7,164,1,1,1,6,2,1,5,3,2,6,1,4,1,1,2,2,1,2,2,7,18,1,2,1,1,2,4,2,4,4,13,4,2,31,2,2,1,2,2,3,6,8,2,2,2,1,2,3,3,3,16,19,1,2,20,3,1,6,2,2,1,45,2,6,1,56,1,1,6,3,1,1,2,2,3,3,2,3,49,28,2,2,74,5,37,25,1,1,20,8,1,8,2,28,1,1,1,1,1,1,9,2,2,62]

#Alpha Feb
#[1,7,26,2,1,1,2,1,1,4,1,2,1,1,19,1,1,7,2,2,4,1,13,3,1,3,2,5,21,1,1,2,1,4,2,5,30,2,1,13,4,3,5,1,2,11,2,3,1,3,3,2,2,1,1,4,1,17,1,8,6,1,3,2,1,2,1,1,1,8,1,1,2,3,5,6,2,1,2,2,29,1,1,2,19,1,1,3,3,2,7,1,8,1,42,1,2,1,1,3,5,1,1,1,5,38,9,7,1,1,20,3,1,3,2,1,2,1,2,2,1,1,140,4,1,1,1,1,1,3,1,1,8,1,4,1,1,2,4,1,1,9,1,1,1,12,1,2,3,1,2,86,5,2,1,16,1,1,4,1,3,37,1,3,1,1,1,1,2,1,5,2,2,3,1,1,6,1,1,1,6,1,2,1,1,1,1,1,2,1,7,1,32,31,2,2,1,1,1,1,1,3,1,3,1,41,9,1,1,2,1,2,8,8,1,65,5,1,24,3,3,1,5,1,61,1,11,1,3,1,4,1,1,1,4,3,8,4,11,5,2,1,1,1,1,3,4,1,3,9,2,9,1,10,1,10,6,2,2,1,1,1,1,12,1,1,2,1,1,2,2,1,1,1,4,1,3,1,2,48,1,1,1,1,2,6,47,1,2,2,1,1,1,1,20,66,4,1,2,1,21,2,3,1,1,2,1,2,1,6,1,1,1,1,1,5,1,1,1,2,4,3,7,1,7,1,4,3,2,1,15,3,2,1,111,12,1,9,2,1,2,4,1,1,1]
#Alpha Nov
#[1, 17, 1, 1, 54, 1, 11, 13, 4, 11, 1, 1, 3, 1, 21, 1, 1, 9, 14, 2, 3, 11, 1, 25, 6, 2, 9, 1, 6, 11, 2, 19, 4, 1, 1, 36, 1, 66, 1, 1, 2, 2, 2, 1, 1, 4, 2, 12, 4, 3, 4, 1, 1, 1, 12, 1, 48, 2, 1, 5, 1, 7, 2, 37, 1, 1, 1, 1, 1, 3, 1, 1, 18, 4, 2, 3]
#Delta Nov 
#[6, 1, 1, 1, 1, 1, 10, 3, 1, 2, 2, 1, 2, 1, 1, 3, 1, 6, 10, 15, 1, 1, 3, 1, 1, 3, 3, 30, 1, 3, 1, 1, 3, 3, 1, 1, 5, 2]
#Europe July
#[5, 1, 1, 3, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 2, 5, 6, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 1, 6, 1, 1, 1, 1, 24, 1, 2, 1, 1, 2, 1, 1, 2, 1, 3, 14, 1, 2, 1, 3, 2, 1]

sizeWithParents = {}
sizesP = []
for k,v in data["Switzerland"].items():
    if 'parent' in v:
        sizeWithParents[k] = len(v['sequences'])
        sizesP.append(len(v['sequences']))
#Alpha Apr 22
#[1,1,1,1,1,1,25,7,1,2,1,1,1,1,2,1,13,4,1,56,1,1,1,5,1,9,10,2,22,3,1,4,1,3,1,2,4,1,3,1,1,3,26,1,1,21,4,3,3,1,13,14,4,1,2,1,1,63,2,2,2,21,1,1,1,5,5,3,1,1,4,1,12,1,5,1,1,14,37,3,43,1,2,8,6,1,6,1,2,1,3,7,3,1,2,1,91,5,1,38,1,15,3,2,4,5,1,1,5,9,12,4,2,1,1,1,5,2,8,2,1,1,1,4,1,1,1,3,1,1,3,4,9,5,21,1,2,8,5,5,1,1,1,54,1,58,5,4,1,16,45,2,2,38,14,9,1,5,2,1,1,4,4,2,1,3,10,18,15,53,113,1,2,1,4,1,10,41,6,49,2,1,13,1,2,1,1,2,2,5,3,1,1,23,19,2,20,14,5,1,1,2,1,1,6,1,1,1,7,2,2,9,1,1,1,9,1,2,6,4,2,1,1,9,6,7,24,5,19,1,5,3,2,7,1,7,1,6,11,4,1,1,1,1,3,1,1,10,7,1,3,6,15,17,1,29,1,8,5,3,2,1,80,1,2,1,4,1,1,1,11,1,10,17,1,1,11,1,14,6,4,2,1,1,4,1,1,2,10,10,2,59,3,1,1,5,15,2,19,22,2,20,1,1,1,1,1,1,1,2,5,3,1,40,35,93,1,1,1,2,2,24,3,2,11,3,2,1,1,2,1,6,1,1,9,2,2,27,1,2,1,3,1,1,1,2,1,2,12,1,2,1,4,1,11,53,2,5,8,1,3,12,13,1,35,1,1,28,5,15,1,1,17,16,1,1,2,10,11,1,1,2,2,18,29,79,2,1,2,2,2,3,2,1,26,1,3,1,1,6,1,2,11,1,2,7,1,4,1,1,1,1,35,1,2,1,1,1,1,134,39,16,12,23,11,1,3,3,177,16,4,1,60,1,22,58,4,18,2,5,5,3,1,255,18,34,19,7,9,1,1,11,3,6,25,2,8,2,1,3,1,8,1,1,11,1,1,8,3,7,12,1,2,1,8,2,15,41,18,8,1,3,6,9,1,1,17,1,1,28,1,46,1,1,3,1,1,12,7,6,2,3,2,8,2,2,21,1,1,4,4,3,1,4,6,1,7,12,1,1,4,1,9,1,14,1,4,2,1,1,19,2,5,17,2,1,3,27,3,12,6,1,4,3,75,5,1,10,1,1,1,7,25,2,2,1,1,2,1,1,1,1,4,1,2,2,1,12,7,84,2,2,61,1,34,1,1,4,1,1,1,2,1,1,1,1,18,1,12,5,3,1,8,7,1,9,1,22,51,1,2,1,3,42,5,1,1,1,2,3,1,1,1,1,1,1,2,1,1,2,23,1,5,1,1,4,11,8,45,54,14,2,3,4,18,22,2,25,2,3,1,4,8,1,2,2,2,2,3,13,1,2,1,7,1]

#Alpha Feb
#[4,2,3,1,1,1,8,1,1,4,12,8,6,4,3,1,1,1,9,23,5,1,5,73,5,2,1,1,1,3,5,4,2,1,3,1,2,3,5,30,1,2,30,2,3,14,11,1,4,1,1,16,1,1,1,1,37,14,1,3,17,2,1,11,1,8,1,5,2,1,5,11,3,3,12,1,1,15,17,4,4,1,20,2,38,5,2,1,10,1,1,39,33,1,2,1,1,2,1,4,7,1,2,4,10,1,1,28,47,24,1,9,42,3,1,2,1,2,16,1,14,1,1,6,19,63,9,4,1,1,3,6,1,17,9,1,2,14,1,10,1,4,1,3,12,1,1,3,2,5,1,6,4,1,1,14,1,4,2,1,3,1,26,39,2,15,18,6,10,1,14,40,3,12,1,1,16,8,1,3,1,9,20,1,8,2,1,1,1,1,1,27,3,2,2,1,4,1,10,1,4,1,1,2,1,1,1,7,1,1,1,4,7,1,4,12,1,24,2,7,1,1,1,59,1,2,1,1,2,2,1,1,5,1,31,2,1,6,1,9,4,1,5,18,2,2,2,38,5,1,9,1,1,2,2,8,55,100,1,1,2,1,1,1,40,8,2,1,20,11,2,2,4,1,10,13,7,8,1,1,2,1,2,1,1,8,1,1,2,1,1,1,1,13,22,8,5,2,7,10,1,4,2,41,1,3,2,1,1,1,9,1,6,3,1,18,1,35,3,1,39,1,1,1,1,45,2,2,2,9,1,7,1,8,29,2,26,1,3,15,1,5,11,1,32,1,3,8,1,2,5,12,53,9,13,1,32,1,1,8,2,1,2,3,8,8,1,2,23,1,1,3,12,1,2,2,2,2,22,2,1,3,2,1,1,2,3,2,2,1,1,16,1,3,1,3,2,1,2,3,1,3,20,1,2,1,8,1,1,1,1,7,12,3,1,3,3,12,7,1,1,2,1,82,1,136,8,4,1,1,4,1,4,2,2,2,3,4,5,9,3,3,1,1,1,1,2,1,27,1,14,5,1,1,1,1,6,5,1,2,1,13,1,1,4,1,1,1,19,7,1,1,1,1,10,1,8,6,1,16,1,11,1,1,1,1,1,1,1,1,21,2,4,1,1,1,2,1,2,1,2,3,1,7,2,1,2,6,4,1,6,17,1,7,1,2,3,2,8,1,2,4,167,39,37,1,6,11,3,24,9,6,9,3,2,3,3,1,3,2,9,4,1,58,1,9,20,2,58,231,18,36,32,7,3,2,1,1,4,4,1,1,5,2,1,2,1,3,3,30,1,3,4,1,3,1,1,2,7,45,1,1,2,29,1,5,36,1,86,1,2,2,1,5,1,1,13,1,2,1,2,3,11,1,1,3,4,6]
#Alpha Nov
#[1, 2, 8, 9, 1, 1, 1, 20, 2, 1, 1, 2, 2, 6, 4, 4, 5, 7, 25, 12, 19, 7, 2, 13, 1, 2, 1, 3, 20, 9, 2, 1, 4, 1, 6, 1, 1, 1, 2, 4, 1, 13, 1, 1, 1, 1, 3, 2, 4, 22, 2, 7, 11, 2, 31, 1, 6, 2, 6, 3, 3, 5, 10, 14, 3, 13, 1, 20, 4, 1, 1, 1, 1, 1, 3, 1, 23, 1, 5, 1, 1, 1, 2, 20, 3, 8, 3, 10, 11, 4, 1, 3, 1, 2, 1, 32, 14, 1, 21, 57, 4, 1, 1, 3, 9, 3, 74, 1, 1, 1, 11, 1, 6, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 10, 1, 8, 2, 5, 2, 8, 2, 1, 25, 40, 110, 25, 14, 4, 3, 2, 135, 33, 20, 2, 1, 1, 1, 2, 1, 3, 1, 2, 2, 11, 3, 5, 4, 4, 1, 1, 1, 12, 28, 2, 1, 29, 2, 1, 1, 6, 12, 3, 2, 3, 41, 1, 1, 1, 2, 2, 2, 3, 16, 4, 4, 34, 3, 1, 1, 3, 2, 3, 2, 5, 1, 3, 1, 26, 1, 1, 1, 1, 5, 3, 4, 2, 10, 2, 4, 23, 25, 27, 3, 29, 28, 4, 2, 25, 1, 1, 5, 6, 1, 1, 15, 18, 2, 4, 46, 4, 7, 9, 1, 1, 2, 7, 1, 8, 4, 1, 4, 1, 3, 2, 23, 2, 1, 4, 1, 1, 1, 1, 89, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 9, 17, 2, 1, 27, 1, 1, 3, 4, 9, 12, 1, 1, 2, 1, 11, 16, 51, 3, 16, 2, 15, 22, 2, 12, 61, 8, 8, 5, 1, 1, 2, 3, 2, 2, 3, 2, 2, 1, 6, 1, 4, 2, 23, 1, 69, 1, 1, 1, 2, 2, 2, 1, 3, 3, 1, 8, 1, 1, 1, 1, 5, 3, 35, 1, 8, 1, 15, 4, 72]
#Delta Nov
#[3, 2, 2, 4, 7, 1, 1, 1, 12, 1, 1, 1, 1, 2, 13, 4, 3, 1, 3, 4, 1, 1, 1, 2, 13, 5, 1, 1, 1, 3, 2, 14, 1, 1, 9, 1, 2, 4, 10, 2, 1, 1, 10, 3, 20, 8, 1, 5, 1, 17, 1, 1, 3, 1, 1, 1, 10, 5, 1, 1, 1, 6, 2, 1, 1, 4, 1, 1, 2, 1, 5, 1, 2, 1, 6, 9, 7, 2, 1, 1, 4, 9]
#Europe July
#[3, 2, 5, 1, 6, 5, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 3, 1, 4, 1, 1, 1, 17, 9, 1, 1, 1, 1, 1, 10, 1, 2, 1, 1, 3, 2, 6, 3, 4]


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


len(parentList)
#Alpha Apr 22 - there are 345 parents for the 655 clusters that have swiss parents
#Delta Apr 22 - there are 379 parents for the 892 clusters that have swiss parents

# Alpah Feb - there are 318 parents for the 608 clusters that have swiss parents
# Alpha Nov - there are 164 parents for the 333 clusters that have swiss parents
# Delta Nov - there are 40 parents for the 82 clusters that have Swiss parents
# Europe July - there are 27 parents for the 39 clusters that have Swiss parents

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
        printNodeAndChildren(0, node, clusterNumber)
        clusterNumber += 1
        print("")


len(rootParents)
len([x for x in rootParents if x in sizeSwissOnly.keys()])
# all of the rootParents are also nodes without Swiss parents (as expected)

#This code below DOESN'T WORK WELL if you have a large Swiss Parent cluster
# Then some non-swiss nodes, then swiss nodes again
# These should be counted as separate introductions but this code will make children
# part of the large swiss parent cluster instead!

#clean out any existing node cluster numbers
#for node in tree.find_clades(order='preorder'):
#    if hasattr(node, "cladeNum"):
#        delattr(node, "cladeNum")
#    if hasattr(node, "clade"):
#        delattr(node, "clade")
#
##number the clusters that have parents:
#clusterNumber = 0
#for par in rootParents:
#    names[par].cladeNum = clusterNumber
#    for child in names[par].find_clades():
#        #don't label *all* children, just those we know have swiss parents
#        #and don't overwrite existing clade numbers!
#        if child.name in sizeWithParents.keys() and not hasattr(child, "cladeNum"):
#            child.cladeNum = clusterNumber
#        if child.name == "NODE_0003128":
#            print(child.cladeNum)
#            print(hasattr(child, "cladeNum"))
#    clusterNumber += 1



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
liberalClusterDF.sort_values(by=["Clade","HeadOfClade"]).to_csv(f'liberalClustersDate-{build}.csv', header=True)
#liberalClusDF_nodate = liberalClusterDF.drop("TotalDatesSwissChildren", axis=1)
#liberalClusDF_nodate.sort_values(by=["Clade","HeadOfClade"]).to_csv(f'liberalClusters-{build}.csv', header=True)


#shortcut to print straight out (unordered)
#(pd.DataFrame.from_dict(data=liberalClusterInfo, orient='index')
#    .to_csv('liberalClusters.csv', header=True))


#### Now re-do for conservative


#Can take the entries with no cluster number as before.
#find the nodes that are heads of cluster and reprocess these, including all children

clusterHeadInfo = liberalClusterDF[liberalClusterDF.notnull().HeadOfClade]
headNodeNumbers = clusterHeadInfo.index.values



conservClusterInfo = defaultdict(dict)

#cycle through all head noes
for k in headNodeNumbers:
    v = data["Switzerland"][k]

    #record clade number
    cl = names[k].cladeNum if hasattr(names[k], "cladeNum") else pd.NA
    #find out if 'top' of clade
    headClusterNode = True if k in rootParents else pd.NA
    if headClusterNode == False:
        print(f"Error! not head of node! {k}")

    directChildren = v['sequences']
    numberDirectChildren = len(directChildren)

    #get all child nodes of this node
    allClusNodes = liberalClusterDF[liberalClusterDF['Clade'] == cl].index.values
    childNodes = [c for c in allClusNodes if c != k]

    clusChilds = []
    for child in childNodes:
        v2 = data["Switzerland"][child]
        clusChilds.extend(v2['sequences'])

    # add all children together
    allChilds = []
    allChilds.extend(directChildren)
    allChilds.extend(clusChilds)
    numberAllChildren = len(allChilds)

    #get information about all children
    datesAllChilds = [datetime.strptime(names[ch].date, "%Y-%m-%d") for ch in allChilds]
    meanDate = get_mean_date(datesAllChilds).strftime("%Y-%m-%d")
    minDate = min(datesAllChilds).strftime("%Y-%m-%d")
    maxDate = max(datesAllChilds).strftime("%Y-%m-%d")
    allSwissDates = [d.strftime("%Y-%m-%d") for d in datesAllChilds]
    #get canton stats
    canton_dict = dict(Counter([names[ch].division for ch in allChilds]))
    child_cantons = ', '.join(['{}: {}'.format(*i) for i in canton_dict.items()])

    #get info about other children
#    all_children = names[k].direct_children 
#    numberAllChildren = len(all_children)
#
#    nonSwissChildren = [c for c in all_children if c.country != "Switzerland"]
#    if nonSwissChildren:
#        countriesNonSwissChildren_dict = dict(Counter([c.country for c in nonSwissChildren]))
#        countriesNonSwissChildren = ', '.join(['{}: {}'.format(*i) for i in countriesNonSwissChildren_dict.items()])
#
#        datesNonSwissChilds = [datetime.strptime(ch.date, "%Y-%m-%d") for ch in nonSwissChildren]
#        meanDateNS = get_mean_date(datesNonSwissChilds).strftime("%Y-%m-%d")
#        minDateNS = min(datesNonSwissChilds).strftime("%Y-%m-%d")
#        maxDateNS = max(datesNonSwissChilds).strftime("%Y-%m-%d")
#
#        #figure out earliest non-swiss Seq
#        [c.country for c in nonSwissChildren if datetime.strptime(c.date, "%Y-%m-%d") == min(datesNonSwissChilds)]
#
#    else:
#        countriesNonSwissChildren = ""
#        meanDateNS = ""
#        minDateNS = ""
#        maxDateNS = ""

    conservClusterInfo[k] = {
        "Clade": cl,
        "HeadOfClade": headClusterNode,
        "NumberSwissChildren": numberAllChildren,
        "TotalNumberChildren": numberAllChildren,

        "MeanDateSwissChildren": meanDate,
        "MinDateSwissChildren": minDate,
        "MaxDateSwissChildren": maxDate,
        "CantonsSwissChildren": child_cantons,

        "NonSwissChildrenCountries": "",
        "MeanDateNonSwissChildren": "",
        "MinDateNonSwissChildren": "",
        "MaxDateNonSwissChildren": "",

        "TotalDatesSwissChildren": allSwissDates,
    }

#convert to dataframe
conservativeClusterDF = pd.DataFrame.from_dict(data=conservClusterInfo, orient='index')

# add those entries that are from before with no (null/NA) cluster number
swissOnlyInfo = liberalClusterDF[liberalClusterDF.isnull().Clade]

# add them together to what we have already
conservClusterDF = pd.concat([swissOnlyInfo, conservativeClusterDF])

#output to file
conservClusterDF.sort_values(by=["Clade","HeadOfClade"]).to_csv(f'conservativeClustersDate-{build}.csv', header=True)




# do alpha & delta
# for more European countries? those with sequencing efforts (germany, swiss, denmark)
# compare timing of imports & numbers

#nonvariant comparison - non-alpha-non-delta introductions

#further incorporate size of cluster? size of swiss-only ness?
#table of list of clusters, dates, size, swiss parent etc - canton

# alpha/delta clusters - get table & builds
# 677
# G5 


#plot a little in R
build = "Alpha"
level = "conservative"
setwd("C:/Users/Emma/wsl/corona/2021_Delta/")
clusDat <- read.csv(paste(level,"Clusters-",build, ".csv",sep=""), as.is=T)

clusDat$MeanDateSwissChildren <- as.Date(clusDat$MeanDateSwissChildren)
clusDat$MaxDateSwissChildren <- as.Date(clusDat$MaxDateSwissChildren)
clusDat$MinDateSwissChildren <- as.Date(clusDat$MinDateSwissChildren)

plot(clusDat$MeanDateSwissChildren, clusDat$NumberSwissChildren,
    ylab="Number of Swiss Sequences", xlab="Date of Sample", main=paste(level, build))
points(clusDat$MinDateSwissChildren, clusDat$NumberSwissChildren, col="red", pch=91)
points(clusDat$MaxDateSwissChildren, clusDat$NumberSwissChildren, col="blue", pch=93)

segments(clusDat$MinDateSwissChildren, clusDat$NumberSwissChildren,
         clusDat$MaxDateSwissChildren, clusDat$NumberSwissChildren, lty="dashed")



