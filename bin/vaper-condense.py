#!/usr/bin/env python3

# vaper-condense.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

version = 1.0

import sourmash
import screed
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from scipy.spatial.distance import squareform
import random
import csv
import random
import datetime
import argparse
import matplotlib.pyplot as plt
import itertools
import re


#----- ARGS -----#
parser = argparse.ArgumentParser()
parser.add_argument("--fasta", dest="fasta", type=str, help="Path to FASTA file.")
parser.add_argument("--stats", dest="read_stats", type=str, help="Path to read stats file.")
parser.add_argument("--dist_threshold", dest="dist_threshold",  default=0.05, type=float, help="Distance thresholds used to create clusters (1-%ANI/100)(default: 0.05)")
parser.add_argument("--ksize", dest="ksize",  default=31, type=int, help="kmer size used by Sourmash")
parser.add_argument("--outdir", dest="outdir",  default='./', type=str, help="Name of outputfile.")
parser.add_argument("--prefix", dest="prefix",  default='', type=str, help="Prefix")
parser.add_argument('--version', action='version', version=f'{version}')
args = parser.parse_args()

start = f"""
vaper-condense.py v{version}

Written by Jared Johnson
jared.johnson@doh.wa.gov

"""
print(start)

#----- FUNCTIONS -----#
# Function to compute the ANI distance for a pair of items
def computeDistance(pair, data):
    id1, id2 = pair
    x = data[id1]['mh'].containment_ani(data[id2]['mh'])
    return { 'id1': id1, 'id2': id2, 'dist': x.dist }

# Function to parallelize the computation of the matrix
def computeMatrix(data):
    dist_long = []

    # Get pairwise comparisons
    ids = list(set([ key for key, value in data.items() ]))
    ids.sort()
    pairs = list(itertools.combinations(ids, 2))
    for pair in pairs:
        dist_long.append(computeDistance(pair, data))

    # write pairwise comparisons to file
    dictToCsv(dist_long, f'{args.outdir}/dists.csv', ["id1","id2","dist"])

    # Extract unique samples 
    ids = list(set(d['id1'] for d in dist_long).union(d['id2'] for d in dist_long)) 
    ids.sort()
    # Initialize the matrix with zeros 
    mat = [[0 for _ in ids] for _ in ids]

    # Fill the matrix with distances from the long format data 
    for entry in dist_long:
        i = ids.index(entry['id1'])
        j = ids.index(entry['id2'])
        mat[i][j] = entry['dist'] 
        mat[j][i] = entry['dist']
    
    mat = np.array(mat)
    np.fill_diagonal(mat, 0)
    mat = squareform(mat)
    
    return mat, ids, dist_long

# Function for assigning clusters using a distance threshold
def createClusters(data, threshold, start, stage):
    result = data
    # Compute pairwise ANI
    mat, ids, dist_long = computeMatrix(data)
   
    # Compute complete-linkage of pairwise ANI
    Z = linkage(mat, method='complete')

    # Cut dendrogram at distance threshold
    max_cluster = start
    clusters = cut_tree(Z, height=threshold).flatten()
    for index, id in enumerate(ids):
        cluster = int(clusters[index]) + int(start) + 1
        result[id]['cluster'] = cluster
        # Calculate the largest cluster number
        if cluster > int(max_cluster):
            max_cluster = cluster
    
    try:
        # Plot the dendrogram
        plt.figure(figsize=(10, 5))
        dendro = dendrogram(Z, labels=ids, orientation='left')

        # Draw a horizontal line at the threshold
        plt.axvline(x=threshold, color='r', linestyle='--')
        plt.title(f'{stage}')
        plt.xlabel('Samples')
        plt.ylabel('Distance')

        # Save the plot to a file
        plt.savefig(f'{args.outdir}/{stage}.png', format='png')
        plt.close()
    except:
        print('Plot not made!')
        
    return result, max_cluster, dist_long
# Function for condensing sequences
def condense(data):
    # group by cluster
    clusters = {}
    for key, value in data.items():
        if value['cluster'] in clusters:
            clusters[value['cluster']].append(key)
        else:
            clusters[value['cluster']] = [ key ]
    # select references for each cluster
    refs = []
    for key, value in clusters.items():
        if len(value) > 1:
            max_covXdepth = 0
            select_ref = []
            for id in value:
                if int(data[id]['covXdepth']) > max_covXdepth:
                    select_ref = [ id ]
                    max_covXdepth = int(data[id]['covXdepth'])
                elif int(data[id]['covXdepth']) == max_covXdepth:
                    select_ref.append(id)
            refs.append(select_ref[0])
        else:
            select_ref = [ value[0] ]
        refs.append(select_ref[0])
        for id in value:
            if id in refs:
                data[id]['condensed'] = ''
            else:
                data[id]['condensed'] = select_ref[0]
    return refs, data

# Function for writing a dictionary to CSV
def dictToCsv(dicts, filename, cols):
    # Write to CSV
    with open(filename, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=cols, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(dicts)

#----- MAIN -----#
# create a dictionary of sequence coverage stats
read_stats = {}
with open(args.read_stats, newline='') as csvfile: 
    csvreader = csv.reader(csvfile, delimiter = '\t')
    header = next(csvreader)
    name_index   = header.index('#rname')
    cov_index    = header.index('coverage')
    depth_index  = header.index('meandepth')
    for row in csvreader:
        name  = row[name_index]
        cov   = row[cov_index]
        depth = row[depth_index]
        read_stats[name] = { 'cov': cov, 'depth': depth }

# Create dictionary of hashes
minhashes = {}
ids = []
for record in screed.open(args.fasta):
    name = re.sub(args.prefix + '_', '', record.name)
    mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=1000)
    mh.add_sequence(record.sequence, True)
    minhashes[ name ] = { 'assembly': record.name, 'mh': mh, 'covXdepth': float(read_stats[ name ]['cov'])/100 * float(read_stats[ name ]['depth']) }
ids = list(set(ids))
print(f'{datetime.datetime.now()}: Loaded {len(minhashes.items())} sequences from {args.fasta}')

# Cluster sequences
print(f'{datetime.datetime.now()}: Creating clusters')
data, max_cluster, dist_long = createClusters(minhashes, args.dist_threshold, 0, 'condense')
print(f'{datetime.datetime.now()}: {len(data.items())} clusters created.')
if len(data.items()) == max_cluster:
    print(f'{datetime.datetime.now()}: Nothing to condense!')
    refs = ids
    for key, value in data.items():
        value['condensed'] = ''
else:
    print(f'{datetime.datetime.now()}: Condensing references.')
    refs, data = condense(data)

# Write results to CSV
print(f'{datetime.datetime.now()}: Writing file.')
results = []
for key, value in data.items():
    del value['mh']
    del value['cluster']
    results.append(value)

dictToCsv(results, f'{args.outdir}/condensed.csv', ["assembly","condensed","covXdepth"])

# End
end = f"""
Results saved to {args.outdir}
"""
print(end)