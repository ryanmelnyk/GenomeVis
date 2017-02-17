#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse, math
import matplotlib.colors as colors
from ete3 import Tree, TreeStyle, NodeStyle
import seaborn as sns

def parse_args():
	parser = argparse.ArgumentParser(description='''
Take a given group from PyParanoid and plot it on a tree, recursively over branches.
	''')
	parser.add_argument('treefile', type=str,help='path to Newick format Tree')
	parser.add_argument('pypdir',type=str,help="path to PyParanoid directory")
	parser.add_argument('group',type=str,help="name of group")
	return parser.parse_args()

def parse_tree(treefile, datadict, group):
	ts = TreeStyle()
	tree = Tree(os.path.abspath(treefile))

	# pal = sns.diverging_palette(10, 240, s=80, l=55, n=13, center="dark")
	pal = sns.cubehelix_palette(rot=-.4, n_colors=13)
	print pal

	for node in tree.iter_descendants("preorder"):

		this_node = []
		nstyle = NodeStyle()
		nstyle["shape"] = "circle"

		if node.is_leaf():
			try:
				if datadict[node.name] > 0:
					nstyle["fgcolor"] = colors.rgb2hex(pal[12])
				else:
					nstyle["fgcolor"] = colors.rgb2hex(pal[0])
			except KeyError:
				nstyle["fgcolor"] = colors.rgb2hex(pal[0])
		else:
			for x in node.iter_descendants("preorder"):
				if x.is_leaf():
					this_node.append(x.name)
			## only generate presence/absence for nodes containing 5+ leaves
			count = 0
			for t in this_node:
				try:
					if datadict[t] > 0:
						count +=1
				except KeyError:
					print t, "not in PyParanoid DB"
			v = int(round(float(count)/float(len(this_node))*12))
			print count, "out of", len(this_node), "\t", v
			nstyle["fgcolor"] = colors.rgb2hex(pal[v])
			nstyle["size"] = 3*math.sqrt(len(this_node))
		node.set_style(nstyle)

	tree.render(group+".png",tree_style=ts,w=12,units="in")
	return

def parse_pyp(pypdir,group):
	for line in open(os.path.join(pypdir,"homolog_matrix.txt"),'r'):
		if line.startswith("\t"):
			header = line.rstrip().split("\t")[1:]
		if not line.startswith(group):
			continue
		else:
			vals = [int(x) for x in line.rstrip().split("\t")[1:]]
	return dict(zip(header,vals))

def main():
	args = parse_args()

	datadict = parse_pyp(os.path.abspath(args.pypdir),args.group)
	# pprint.pprint(datadict)
	parse_tree(args.treefile, datadict, args.group)


if __name__ == '__main__':
	main()
