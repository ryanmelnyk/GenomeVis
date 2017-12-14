#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

from ete3 import Tree, TreeStyle, PieChartFace,NodeStyle
import seaborn as sns
import os
import numpy as np
from Bio import SeqIO
import matplotlib.colors as colors
import argparse

def parse_args():
	parser = argparse.ArgumentParser(description='''
Takes a tree and generates an image to show the general degree of core genome
conservation among different clades.
	''')
	parser.add_argument('tree', type=str,help='path to the tree file')
	parser.add_argument('pypdir',type=str,help="path to PyParanoid directory")
	parser.add_argument('genomedb',type=str,help="path to genome database folder")
	return parser.parse_args()

def add_piecharts(node,nstyle,strains,a,genome_sizes,pal):
	node_strains = node.get_leaf_names()
	node.set_style(nstyle)
	if len(node_strains) > 1:
		col_indices = [strains.index(ns) for ns in node_strains]
		orth_count = 0
		for i in range(a.shape[0]):
			if np.count_nonzero(a[i,col_indices]) == len(col_indices):
				orth_count += 1
		frac = float(orth_count)/float(sum([genome_sizes[ns] for ns in node_strains])/len(node_strains))

		size = 15
		f = PieChartFace([frac*100,100-(frac*100)],size,size,colors=[colors.rgb2hex(c) for c in [pal[0],pal[7]]])
		node.add_face(f,0)
	return

def main():
	args = parse_args()
	tree = Tree(os.path.abspath(args.tree))
	strains = tree.get_leaf_names()

	matrixfile = open(os.path.join(os.path.abspath(args.pypdir),"homolog_matrix.txt"),'r')
	header = matrixfile.readline().rstrip().split("\t")
	indices = [header.index(s) for s in strains]

	lines = []
	for line in matrixfile:
		vals = line.rstrip().split("\t")
		lines.append([int(bool(int(vals[i]))) for i in indices])
	a = np.stack(lines)

	print a.shape

	genome_sizes = {}
	for s in strains:
		count = 0
		for seq in SeqIO.parse(open(os.path.join(os.path.abspath(args.genomedb),"pep/{}.pep.fa".format(s)),'r'),'fasta'):
			count += 1
		genome_sizes[s] = count

	pal = sns.color_palette("Set2",10)
	sns.palplot(pal)
	nstyle = NodeStyle()
	nstyle["size"] = 0
	ts = TreeStyle()
	for node in tree.iter_descendants("preorder"):
	    add_piecharts(node,nstyle,strains,a,genome_sizes,pal)
	add_piecharts(tree,nstyle,strains,a,genome_sizes,pal)

	namedict = {l[0] : l[1] for l in [line.rstrip().split("\t") for line in open("renamed_list.txt",'r').readlines()]}

	for node in tree.iter_descendants("preorder"):
		if node.is_leaf():
			node.name = namedict[node.name]
	tree.render("test.pdf",tree_style=ts,h=96,units="in")


if __name__ == '__main__':
	main()
