#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Given a reference genome file in genbank format, use homolog matrix information
to map regions found in other strains. Should be very helpful in identifying
genomic islands in otherwise closely related strains.
	''')
	parser.add_argument("gbk",type=str,help="genbank file containing locus order")
	parser.add_argument('map_to', type=str,help='strain to map to')
	parser.add_argument('strains',type=str,help="list of strains to check for presence/absence")
	parser.add_argument('pyp',type=str,help="directory containing PyParanoid output")
	parser.add_argument('outfile',type=str,help="path to write output to")
	return parser.parse_args()

def get_locus_tags(gbk):
	seqs = {}
	for seq in SeqIO.parse(open(gbk,'r'),'genbank'):
		seqs[seq.id] = []
		for feat in seq.features:
			if feat.type == "CDS":
				try:
					seqs[seq.id].append(feat.qualifiers["locus_tag"][0])
				except KeyError:
					seqs[seq.id].append(feat.qualifiers["protein_id"][0])
	return seqs

def get_groups(pyp, seqs, map_to):
	groupdict = {}

	locustags = []
	for seq in seqs:
		for s in seqs[seq]:
			locustags.append(s)

	count = 0

	header = open(os.path.join(pyp,"locustag_matrix.txt"),'r').readline().rstrip().split("\t")
	i = header.index(map_to)
	for line in open(os.path.join(pyp,"locustag_matrix.txt"),'r'):
		vals = [v.split(".")[0] for v in line.rstrip().split("\t")[i].split(";")]
		for s in locustags:
			if s in vals:
				groupdict[s] = line.rstrip().split("\t")[0]
		if count % 1000 == 0:
			print count, "lines parsed..."
		count += 1

	return groupdict

def find_homologs(seqs, groupdict, strains, outfile, pyp, map_to):
	header = open(os.path.join(pyp,"locustag_matrix.txt"),'r').readline().rstrip().split("\t")
	indices = [header.index(x) for x in strains]
	groups = [g for g in groupdict.values()]
	datadict = {}
	for line in open(os.path.join(pyp,"locustag_matrix.txt"),'r'):
		vals = line.rstrip().split("\t")
		if vals[0] in groups:
			datadict[vals[0]] = []
			for i in indices:
				if len(vals[i].split(";")) > 1:
					datadict[vals[0]].append("Multiple")
				else:
					datadict[vals[0]].append(vals[i])

	group_annotations = get_group_annotations(pyp)

	o = open(outfile,'w')
	for seq in seqs:
		o.write(">{}\n".format(seq))
		o.write("{}\tgroup\t{}\tannotation\n".format(map_to,"\t".join(strains)))
		for s in seqs[seq]:
			try:
				o.write("{}\t{}\t{}\t{}\n".format(s,groupdict[s],"\t".join(datadict[groupdict[s]]),group_annotations[groupdict[s]]))
			except KeyError:
				pass
	return

def get_group_annotations(pyp):
	group_annotations = {}
	for line in open(os.path.join(pyp,"group_descriptions.txt"),'r'):
		vals = line.rstrip().split("\t")
		counts = {}
		for x in set(vals[1:]):
			counts[x] = vals.count(x)
		group_annotations[vals[0]] = sorted(counts.items(), reverse=True, key = lambda x: x[1])[0][0]
	return group_annotations

def main():
	args = parse_args()
	gbk = args.gbk
	map_to = args.map_to
	strains = [x.rstrip() for x in open(os.path.abspath(args.strains),'r')]
	pyp = os.path.abspath(args.pyp)
	outfile = os.path.abspath(args.outfile)

	seqs = get_locus_tags(gbk)
	groupdict = get_groups(pyp, seqs, map_to)
	find_homologs(seqs, groupdict, strains, outfile, pyp, map_to)

if __name__ == '__main__':
	main()
