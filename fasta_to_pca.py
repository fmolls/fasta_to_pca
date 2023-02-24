import os
from Bio.SeqIO import parse
from sklearn.decomposition import PCA
from itertools import product
import sys
import matplotlib.pyplot as plt
import argparse
from sklearn.preprocessing import StandardScaler
xco = []
yco = []
iterator = 0
fastaData = {}
matrixpca = [[]]
fasta_names = []
pm_pca = {}
parser = argparse.ArgumentParser("This programm calculates the k-mers of choosen Fasta Data and plots a 2d pca. FASTA Data has to be saved to fasta directory and number of k-mers has to be given as command argument -k.")
parser.add_argument("-k","--k_mers", help="number of k-mers to be considered", type=int, choices=[1,2,3,4,5,6,7,8,9])
parser.add_argument("-f","--fasta", nargs='+', help="fasta data for pca")
parser.add_argument("-o","--output", help="if set, plot will be saved to png with parsed name, else plot will be shown")
parser.add_argument("-p", "--position", help='Additional parameter for position marking, if set a additional pca will be generated only with seqences containing set Marker. The first Sequence Position is 0, not 1. This means the Median of a Seqence of length 9 is 4, not 5.')
parser.add_argument("-n","--seqineq", help="If the length of sequences of one Fasta-data are unequal, normalisation will be calculated differently. 0 = equal length of seqeunces , 1 = unequal length of sequences. If 1 is true, Runtime will increase.")
args = parser.parse_args()
def kmer_dict_calc(k):
	# anlegen eines dictionaries für jedes mögliche k-mere, jeweils initialisiert mit 0 -- 1 O(1)
	kmer_dict = {}
	a = "ATGC"
	kmer_perm = product(a, repeat=k)
	for i in list(kmer_perm):
		s = ""
		for l in i:
			s = s + l
		kmer_dict[s] = 0
	return kmer_dict
# -- 1
def normalisation(k, iterator):
	# berechnung des normalisierung O(1) / O(n)
	for filename in args.fasta:
		fastaData[filename] = parse(filename, 'fasta')
	Valuen = list(fastaData[fasta_names[iterator]])
	if int(args.seqineq)==0:
		return len(Valuen) * (len(Valuen[0].seq) - k + 1)
	else:
		intermediate = 0
		for record in Valuen:
			intermediate = len(record.seq) + intermediate
	return len(Valuen) * ((intermediate/len(Valuen)) - k + 1)
def k_mer_calc(k,iterator,kmer_dict):
	# O(n)
	result = {}
	i = 0
	# Algorithmus um die k-mere aus den fastadateien auszulesen und deren anzahl im dictionary zu speichern -- 2
	for Value in fastaData.values():
		kmer_dict_intermediate_pm = kmer_dict.copy()
		kmer_dict_intermediate = kmer_dict.copy()
		norm = 1/normalisation(k, iterator)
		iterator = iterator + 1
		try:
			p = int(args.position)
		except TypeError:
			print("no position marker set")
			# dies iteriert über den einzelnen record und speichert die anzahl der k-mere ab -- 3
		for record in Value:
			start = 0
			end = start + k
			while start < (len(record.seq)-k+1):
				if "N" not in record.seq[start:end]:
					kmer_dict_intermediate[record.seq[start:end]] = kmer_dict_intermediate.get(record.seq[start:end]) + norm
					start = start + 1
					end = end + 1
				else:
					start = start + 1
					end = end + 1
			# -- 3
			# das selbe wie 3, nur das des die k-mere welche auf dem positionsmarker liegen separat abspeichert -- 4
			try:
					if p >= k:
						start_pm = p - k
						end_pm = p
					else:
						start_pm = 0
						end_pm = k
					while start_pm < (len(record.seq) - k + 1) and start_pm < p:
						if "N" not in record.seq[start:end]:
							kmer_dict_intermediate_pm[record.seq[start_pm:end_pm]] = kmer_dict_intermediate_pm.get(
								record.seq[start_pm:end_pm]) + norm
							start_pm = start_pm + 1
							end_pm = end_pm + 1
						else:
							start_pm = start_pm + 1
							end_pm = end_pm + 1
			except UnboundLocalError:
				pass
			# -- 4
		result[fasta_names[i]] = kmer_dict_intermediate
		pm_pca[fasta_names[i]] = kmer_dict_intermediate_pm
		i = i + 1
	return result
	# -- 2
def plot_pca(matrixpca):
	# O(n)
	intermediate = []
	i = 0
	for l in list(matrixpca):
		intermediate.append(list(matrixpca[fasta_names[i]].values()))
		i = i + 1
	scaler = StandardScaler()
	scaler.fit(intermediate)
	scaled_data = scaler.transform(intermediate)
	pca = PCA(n_components=2)
	pca.fit(scaled_data)
	x_pca = pca.transform(scaled_data)
	scaled_data.shape
	x_pca.shape
	plt.scatter(x_pca[:,0],x_pca[:,1])
	for i, txt in enumerate(fasta_names):
		plt.annotate(os.path.basename(txt), (x_pca[i,0],x_pca[i,1]))
	try:
		plt.savefig(args.output)
	except NameError:
		plt.show()
	try:
		intermediate = []
		i = 0
		for l in list(pm_pca):
			intermediate.append(list(pm_pca[fasta_names[i]].values()))
			i = i + 1
		scaler.fit(intermediate)
		scaled_data = scaler.transform(intermediate)
		pca = PCA(n_components=2)
		pca.fit(scaled_data)
		x_pca = pca.transform(scaled_data)
		scaled_data.shape
		x_pca.shape
		plt.scatter(x_pca[:, 0], x_pca[:, 1])
		i = 0
		for i, txt in enumerate(fasta_names):
			plt.annotate(os.path.basename(txt), (x_pca[i,0], x_pca[i,1]))
		try:
			plt.savefig(args.output + "_with_position_marker=" + args.position)
			sys.exit()
		except NameError:
			plt.show()
		except TypeError:
			print("position marker position marker not set")
	except UnboundLocalError:
		print("position marker position marker not set")
		sys.exit()
	else:
		sys.exit()
def main():
	try:
		k = int(args.k_mers)
		kmer_dict = kmer_dict_calc(k)
	except TypeError:
		print("please input number of k-meres ")
		sys.exit()
	else:
			for filename in args.fasta:
				fastaData[filename] = parse(filename,'fasta')
				fasta_names.append(filename)
			plot_pca(k_mer_calc(k,iterator,kmer_dict))


if __name__ == "__main__":
	#O(1) + O(1)/O(n) + O(n) + O(1) = O(n)
    main()	
 


