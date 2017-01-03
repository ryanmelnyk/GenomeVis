#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import random
import argparse, os, sys
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
import seaborn as sns
from Bio.SeqFeature import FeatureLocation

def parse_args():
	parser = argparse.ArgumentParser(description='''
Here we want to plot the regions around 2 homologous genes in different organisms
and color-code surrounding genes based on homology.
	''')
	parser.add_argument('infile',type=str, help='input file containing path to genbank and name of locus and genbank file type')
	parser.add_argument('span',type=int,help='distance (in bp) to draw around central gene')
	parser.add_argument('name',type=str, help='name of plot for output file')
	parser.add_argument('locus_mat',type=str,help="path to matrix file containing locus tag information")
	parser.add_argument('--no_labels',action='store_true',help='use to hide labels')
	return parser.parse_args()

def parse_genbank(g):
	if g[2] == "ensembl":
		FIELD = "protein_id"
	elif g[2] == "prokka":
		FIELD = "locus_tag"
	for seq in SeqIO.parse(open(g[0],'r'),"genbank"):
		for feat in seq.features:
			if feat.type == "CDS":
				if feat.qualifiers[FIELD][0] == g[1]:
					print "Found", g[1], "in", seq.id
					print "loc:",feat.location
					return seq, (int(feat.location.start), int(feat.location.end))

	print "locus not found. try again."
	sys.exit()
	return

def plot(seq, span, coords, g, GD, count, locus_tags, labels):

	if g[2] == "ensembl":
		FIELD = "protein_id"
	elif g[2] == "prokka":
		FIELD = "locus_tag"

	track = GD.new_track(count, height=1, name="CDS",\
		scale_ticks=False,scale_largeticks=False,scale_largetick_labels=False, scale_largetick_interval=10000,\
		scale_smallticks=False, scale_smalltick_labels=False, scale_smalltick_interval=1000,\
		greytrack=False, hide=False, scale=False
	)
	feature_set = track.new_set()


	count = 0
	for feat in seq.features:
		if feat.type == "CDS":
			if int(feat.location.start) > (coords[0]-(span/2)) and int(feat.location.end) < (coords[1]+(span/2)):
				newloc = FeatureLocation(feat.location.start-(coords[0]-(span/2)),feat.location.end-(coords[0]-(span/2)),strand=feat.strand)
				feat.location = newloc
				print newloc
				if g[2] == "prokka":
					feature_set.add_feature(feat, sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=.4,color="#D3D3D3", \
						label=labels,name=feat.qualifiers['product'][0],label_strand=1,label_size = 8,label_position="middle", label_angle=20, \
						border=colors.black)
				feature_set.add_feature(feat, sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=.4,color="#D3D3D3", \
					label=True,name=feat.qualifiers[FIELD][0],label_strand=-1,label_size = 8,label_position="middle", label_angle=90, \
					border=colors.black)
				locus_tags[g[0].split(".")[0]].append(feat.qualifiers[FIELD][0])
		elif feat.type == "gene":
			if g[2] == "ensembl":
				if int(feat.location.start) > (coords[0]-(span/2)) and int(feat.location.end) < (coords[1]+(span/2)):
					newloc = FeatureLocation(feat.location.start-(coords[0]-(span/2)),feat.location.end-(coords[0]-(span/2)),strand=feat.strand)
					feat.location = newloc
					try:
						feature_set.add_feature(feat, sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=.4,color="#D3D3D3", \
							label=labels,name=feat.qualifiers['note'][0],label_strand=1,label_size = 8,label_position="middle", label_angle=20, \
							border=colors.black)
					except KeyError:
						pass

	return

def write(GD,name,span):
	GD.draw(x=0.1,format="linear",\
			orientation = "landscape", \
			track_size=0.035,\
			fragments=1,
			start=0, end=span+2000
	)
	GD.write(name, "PDF")
	# GD.write(name+".svg", "SVG")
	return

def find_homologs(GD, locus_tags, locus_mat):
	groups = {}
	count = 0
	for line in open(locus_mat,'r'):
		vals = [v.split(".")[0] for v in [x.split(";") for x in line.rstrip().split("\t")] for v in v]
		found = []
		for strain in locus_tags:
			for t in locus_tags[strain]:
				if t in vals:
					found.append((strain,t))

		if len(found) > 1:
			for f in found:
				groups[f[1]] = count
			count += 1
	print groups
	return groups

def change_colors(GD, groups):
	cl = [colors.HexColor(c) for c in sns.cubehelix_palette(len(set(groups.values())),dark=0.1,light=0.9,rot=2.5).as_hex()]
	random.shuffle(cl)
	for t in GD.get_tracks():
		for s in t.get_sets():
			for feat in s.get_features():
				if feat.name in groups:
					feat.color = cl[groups[feat.name]]
	return

def main():
	args = parse_args()
	genfiles = [line.rstrip().split("\t") for line in open(os.path.abspath(args.infile),'r')]
	span = args.span
	name = args.name
	locus_mat = os.path.abspath(args.locus_mat)
	if args.no_labels:
		labels = False
	else:
		labels = True

	GD = GenomeDiagram.Diagram('gbk',name)
	count = 1
	locus_tags = {}
	for g in genfiles:
		if g[0].split(".")[0] not in locus_tags:
			locus_tags[g[0].split(".")[0]] = []
		contigseq, coords = parse_genbank(g)
		plot(contigseq, span, coords, g, GD, count, locus_tags, labels)
		count += 1
	groups = find_homologs(GD, locus_tags, locus_mat)
	change_colors(GD, groups)
	write(GD,name,span)


if __name__ == '__main__':
	main()
