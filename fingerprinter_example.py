from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np
import pandas as pd
import matplotlib

smiles = pd.read_csv('smiles.csv', sep=',', header=None)

compounds = []
fps = []

#convert smiles to rkd objects
for i in range(smiles.shape[0]):	
	compounds.append(Chem.MolFromSmiles(smiles.ix[i, 0]))

#draw all structures to png file
img = Draw.MolsToGridImage(compounds, molsPerRow=4, subImgSize=(200, 200))
img.save('structures.png')
	
#compute fingerprinters
for m in compounds:
	arr = np.zeros((1, ))
	fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024)
	DataStructs.ConvertToNumpyArray(fp, arr)
	fps.append(arr)

#create dataframe to store fingerprinters and write to csv file
df = pd.DataFrame(fps, index=smiles.ix[:, 0])
df.to_csv('fingerprinters.csv')


