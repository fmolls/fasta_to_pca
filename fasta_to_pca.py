import os
from Bio.SeqIO import parse
from Bio.Cluster import pca
from itertools import product
import sys
import matplotlib.pyplot as plt
import argparse
xco = []
yco = []
fastaData = {}
matrixpca = [[]]
fasta_names = []
pm_pca = {}
parser = argparse.ArgumentParser("This programm calculates the k-mers of choosen Fasta Data and plots a 2d pca. FASTA Data has to be saved to fasta directory and number of k-mers has to be given as command argument -k.")
parser.add_argument("-k","--k_mers", help="number of k-mers", type=int, choices=[1,2,3,4,5,6,7,8,9])
parser.add_argument("-f","--fasta", nargs='+', help="fasta data for pca")
parser.add_argument("-o","--output", help="if set, plot will be saved to png with parsed name, else plot will be")
parser.add_argument("-p", "--position", help='Additional parameter for position marking, if set a additional pca will be generated only with seqences containing set Marker')
args = parser.parse_args()
    
def k_mer_calc(k):
	kmer_dict = {}
	a = "ATGC"
	kmer_perm = product(a, repeat = k)
	for i in list(kmer_perm):
		s = ""
		for l in i:
			s = s + l
		kmer_dict[s] = 0
	result = {}
	i = 0
	for Value in fastaData.values():
		kmer_dict_intermediate_pm = kmer_dict.copy()
		kmer_dict_intermediate = kmer_dict.copy()
		try:
			p = int(args.position)
		except TypeError:
			print("no position marker set")
		for record in Value:
			start = 0
			end = start + k
			while start < (len(record.seq)-k+1):
				kmer_dict_intermediate[record.seq[start:end]] = kmer_dict_intermediate.get(record.seq[start:end]) + 1
				start = start + 1
				end = end + 1
			try:
					if p >= k:
						start_pm = p - k
						end_pm = p
					else:
						start_pm = 0
						end_pm = k
					while start_pm < (len(record.seq) - k + 1) and start_pm < p:
						kmer_dict_intermediate_pm[record.seq[start_pm:end_pm]] = kmer_dict_intermediate_pm.get(
							record.seq[start_pm:end_pm]) + 1
						start_pm = start_pm + 1
						end_pm = end_pm + 1
			except UnboundLocalError:
				pass
		result[fasta_names[i]] = kmer_dict_intermediate
		pm_pca[fasta_names[i]] = kmer_dict_intermediate_pm
		i = i + 1
	return result

def plot_pca(matrixpca):
	intermediate = []
	i = 0
	for l in list(matrixpca):
		intermediate.append(list(matrixpca[fasta_names[i]].values()))
		i = i + 1
	columnmean, coordinates, components, eigenvalues = pca(intermediate)
	i = 0
	for filename in args.fasta:
		xco.append(coordinates[i][0])
		yco.append(coordinates[i][1])
		i = i + 1
	plt.scatter(xco, yco)
	i = 0
	for i, txt in enumerate(fasta_names):
		plt.annotate(txt, (xco[i],yco[i]))
	try:
		plt.savefig(args.output)
	except NameError:
		plt.show()
	try:
		p = int(args.position)
		xco_pm = []
		yco_pm = []
		intermediate = []
		i = 0
		for l in list(pm_pca):
			intermediate.append(list(pm_pca[fasta_names[i]].values()))
			i = i + 1
		columnmean_pm, coordinates_pm, components_pm, eigenvalues_pm = pca(intermediate)
		i = 0
		for filename in args.fasta:
			xco_pm.append(coordinates_pm[i][0])
			yco_pm.append(coordinates_pm[i][1])
			i = i + 1
		plt.scatter(xco_pm, yco_pm)
		i = 0
		for i, txt in enumerate(fasta_names):
			plt.annotate(txt, (xco_pm[i], yco_pm[i]))
		try:
			plt.savefig(args.output + "_with_position_marker=" + args.position)
			sys.exit()
		except NameError:
			plt.show()
	except UnboundLocalError:
		print("position marker position marker not set")
		sys.exit()
	else:
		sys.exit()

    
def main():
	try:
		k = int(args.k_mers)
	except TypeError:
		print("please input number of k-meres ")
		sys.exit()
	else:
		try:
			for filename in args.fasta:
				fastaData[filename] = parse(os.getcwd() + "/" + filename,'fasta')
				fasta_names.append(filename)
			plot_pca(k_mer_calc(k))
		except TypeError:
			pass
			sys.exit()
		else:
			sys.exit()	

if __name__ == "__main__":
    main()	
 


