#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
import seaborn as sns
from ete3 import Tree
from itertools import combinations
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

def parse_args():
	parser = argparse.ArgumentParser(description='''
Using seaborn to plot heatmap of homologs for a given tree.
	''')
	parser.add_argument('tree_loc', type=str,help='path to tree file')
	parser.add_argument('outdir',type=str,help='PyParanoid output directory')
	parser.add_argument('groups',type=str,help='file containing groups and their annotation to display')
	parser.add_argument('--dpi',type=int,help="dpi")
	return parser.parse_args()

def dist_from_tree(tree_loc):
	tree = Tree(tree_loc)
	leaves = tree.get_leaf_names()
	print leaves
	return leaves

def make_matrix(outdir,strains,groups):
	dat = {}
	for g in groups:
		present = []
		for seq in SeqIO.parse(open(os.path.join(outdir,"homolog_faa",g[0]+".faa"),'r'),'fasta'):
			if str(seq.id).split("|")[0] not in present:
				present.append(str(seq.id).split("|")[0])
		for seq in SeqIO.parse(open(os.path.join(outdir,"prop_homolog_faa",g[0]+".faa"),'r'),'fasta'):
			if str(seq.id).split("|")[0] not in present:
				present.append(str(seq.id).split("|")[0])
		dat[g[1]] = pd.Series(dict(zip(strains,[present.count(s) for s in strains])))
	df = pd.DataFrame(dat).reindex(strains)
	cols = [x[1] for x in groups]
	return df[cols]

def plot_clustergrid(df,prefix,groups,dpi):
	fig, ax = plt.subplots()
	hm = sns.heatmap(df, linewidths=.4,cbar=False,cmap="Greens",ax=ax,square=True,xticklabels=["" for x in range(df.shape[1])],yticklabels=["" for x in range(df.shape[0])])
	ax.yaxis.tick_right()
	plt.yticks(rotation=0)
	plt.xticks(rotation=-40)
	plt.savefig('{}.png'.format(prefix),format='png',dpi=400)
	plt.savefig('{}.svg'.format(prefix),format='svg')
	return

def main():
	args = parse_args()
	if args.dpi:
		dpi = args.dpi
	else:
		dpi = 300
	tree_loc = os.path.abspath(args.tree_loc)
	outdir = os.path.abspath(args.outdir)
	groups = [(v[0],v[1]) for v in [line.rstrip().split("\t") for line in open(os.path.abspath(args.groups),'r')]]
	print groups
	leaves = dist_from_tree(tree_loc)
	df = make_matrix(outdir,leaves,groups)
	plot_clustergrid(df,os.path.basename(args.groups).split(".")[0],groups,dpi)



if __name__ == '__main__':
	main()
