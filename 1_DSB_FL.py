# first attempt - fingerprint similarity comparison, no 3D structure 
import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
import pandas as pd

# convert SMI to csv - use pandas to convert text file to dataframe 
df = pd.read_csv('/home/martyn0000/fragrance_investigation/GDB-13-FL.smi',
                sep='\r\n', engine='python')
df = df.iloc[1: , :]
print(df.head(2))

# sanitise input - smiles are canonical, do they need aromatised? 

# rdkit calculate features - Morgan fingerprint 
# add column heading 'Smiles' 
PandasTools.AddMoleculeColumnToFrame(df,'canonical_smiles','molecule',includeFingerprints=True)
df['morgan'] = df['molecule'].map(lambda x:GetMorganFingerprintAsBitVect(x,2)) 

# similarity search of database 
similarity_target=Chem.MolFromSmiles("") 
def get_dataframe_with_x_most_similar_compounds_to_query(query, mol_df, molCol='molecule', x=20):
    query_fp = GetMorganFingerprintAsBitVect(query,2)
    mol_df['similarity'] = mol_df['morgan'].map(lambda x:DataStructs.TanimotoSimilarity(query_fp, x))
    mol_df.sort_values(['similarity'], ascending=False, inplace=True)
    return mol_df[:x].copy()  
similarity_df_filter=get_dataframe_with_x_most_similar_compounds_to_query(similarity_target, df)
similarity_df_filter.head()

# generate similarity maps of most suitable candidates 

# file abandoned - need to remake as a jupyter notebook due to the size of the dataframe, can't repeatedly rerun all the code 