import pandas as pd
import os
from data_preprocessing import preprocess_dataframe
from descriptor_calculation import descriptor_cal

PATH= "All Properties with Finguprints_3.csv"
descriptor_list= ["Morgan fingerprints","MACCSKeysFingerprint","RDKitDescriptors",'Mordred descriptors','PubChemFingerprint'] 
# descriptor_list= ["MACCSKeysFingerprint","RDKitDescriptors"] 
df_size= 10

# read data
df = pd.read_csv(PATH)
# drop duplicates and encodes taget
df = preprocess_dataframe(df) # smiles , target 
print(f"shape of df after preprocessing: {df.shape}")
# Loop through the list of descriptors
for descriptor in descriptor_list:
    # Call the descriptor_cal function to get X and Y
    X, Y = descriptor_cal(df.head(df_size), descriptor)
    
    # Create a folder for the descriptor if it doesn't exist
    folder_path = os.path.join(os.getcwd(),"Descriptors", f"{descriptor}")
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    # Save the X and Y DataFrames to separate CSV files in the folder
    x_file_path = os.path.join(folder_path, 'X.csv')
    y_file_path = os.path.join(folder_path, 'Y.csv')
    X.to_csv(x_file_path, index=False)
    Y.to_csv(y_file_path, index=False)

    print(f"X and y saved for{descriptor}")
    print(f"descriptors shape: {X.shape}")



