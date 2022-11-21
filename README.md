# fasta_to_pca
 A Tool for visualizing the k-mers of different fasta-datas as PCA. Can be used as quality-control of samples.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

The source code can be cloned to your local directory using:
```
git clone https://github.com/fmolls/bioinf_projekt_malewins.git
```

or downloaded and extracted from the github page: [https://github.com/malewins/SEQing](https://github.com/fmolls/bioinf_projekt_malewins)


The file ```requirements.txt``` can be used to install all necessary dependencies. Python 3.8 or higher is required and we recommend to setup a virtual environment for this project. If your current python points to a python2 version, please put ```python3``` instead of just ```python``` before running SEQing. The same applies to the package installer ```pip```.

Once you have setup your virtual environment run the following code to install the dependencies:
```
pip install -r requirements.txt
```
or 
```
pip3 install -r requirements.txt
```
if your pip points to an existing Python2 environment.

## Running fasta_to_pca with sample data

To test-run fasta_to_pca with sample data, you have to create a directory named fastas in your working directory and place the example.fasta in it.
Then you test-run the programm:

```
python3 seq.py -k 3
```
the result should look like 

![example](example.png)

## Running fasta_to_pca

To run fasta_to_pca, you have to create a fastas directory in your working directory. Now you have to place every fasta-data you want to compare in this directory. Now you have to choose the number of k-mers and run this command in your working directory:

```
python3 seq.py -k 3
```

In this example we choose the k-mer number 3. After compiling a png with the pca of your fastas should pop up.

