from tkinter import *
from tkinter import ttk
from PIL import Image, ImageTk

from rdkit import Chem
from rdkit.Chem import Draw

import pandas as pd
idx = 0
compound_name = ""

def int_with_default(value, default):
    res = 0
    try:
        res = int(value)
    except ValueError:
        res = default
    return res

def search_compound():
    compound_row = df.loc[df['Smiles']==txt_compound.get()].to_numpy()[0]
    compound = compound_row[0]
    if (compound == ''):
        compound = txt_compound.get()

    img = Draw.MolToImage(Chem.MolFromSmiles(compound))
    img = ImageTk.PhotoImage(img)

    lbl_molecule.configure(image=img)
    lbl_molecule.image=img

    absorption_max = int_with_default(compound_row[2], -1)
    emission_max = int_with_default(compound_row[3], -1)
    if (absorption_max == -1 or emission_max == -1):
        lbl_status.configure(text="NOT Fluroscent")
    else:
        lbl_status.configure(text="Fluroscent")
