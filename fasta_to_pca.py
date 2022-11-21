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
parser = argparse.ArgumentParser("This programm calculates the k-mers of choosen Fasta Data and plots a 2d pca. FASTA Data has to be saved to fasta directory and number of k-mers has to be given as command argument -k.")
parser.add_argument("-k","--k_mers", help="number of k-mers", type=int, choices=[1,2,3,4,5,6,7,8,9])
parser.add_argument("-f","--fasta", nargs='+', help="fasta data for pca")
parser.add_argument("-o","--output", help="if set, plot will be saved to png with parsed name, else plot will be")
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
		kmer_dict_intermediate = kmer_dict.copy()
		for record in Value:
			start = 0
			end = start + k
			while start < (len(record.seq)-k+1):
				kmer_dict_intermediate[record.seq[start:end]] = kmer_dict_intermediate.get(record.seq[start:end]) + 1
				start = start + 1
				end = end + 1
		result[fasta_names[i]] = kmer_dict_intermediate
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
		sys.exit()
	except NameError:
		plt.show()
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
			print("please input fasta data for computing pca")
			sys.exit()
		else:
			sys.exit()	

if __name__ == "__main__":
    main()	
 


