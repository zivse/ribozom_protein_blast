import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px


def get_protein_by_place(place, human_refrence):
    values = human_refrence.values()
    for index, value in enumerate(values):
        if place < value:
            return index - 1


rows = 75
cols = 75

# Create a 75x75 matrix filled with zeros
counting_matrix = [[0 for _ in range(cols)] for _ in range(rows)]


human_refrence = {
    "XXXXXXXXXXXXX": 0,
    "Q9BZE1": 423,
    "Q9BRJ2": 729,
    "Q96DV4": 1109,
    "Q8TCC3": 1270,
    "O60783": 1398,
    "Q9H0U6": 1578,
    "P82930": 1796,
    "Q9P015": 2092,
    "Q9NZE8": 2280,
    "Q9Y2Q9": 2467,
    "P82909": 2570,
    "Q9Y3D9": 2760,
    "O15235": 2898,
    "O95886": 3877,
    "P82921": 3964,
    "Q9BYC8": 4152,
    "Q9BYD3": 4463,
    "P82914": 4720,
    "Q9NWU5": 4926,
    "Q6P161": 5064,
    "Q96A35": 5280,
    "Q9NYK5": 5618,
    "Q9Y291": 5724,
    "Q96EL3": 5836,
    "P82675": 6266,
    "Q92665": 6661,
    "P52815": 6859,
    "Q13084": 7115,
    "Q9H9J2": 7447,
    "Q8N983": 7662,
    "Q7Z2W9": 7867,
    "P82932": 7992,
    "Q9BYD1": 8170,
    "P82933": 8566,
    "Q9Y6G3": 8708,
    "Q16540": 8861,
    "Q8IXM3": 8998,
    "Q9Y3B7": 9190,
    "Q9NQ50": 9396,
    "Q9Y2R9": 9638,
    "P09001": 9986,
    "Q9P0M9": 10134,
    "Q9Y399": 10430,
    "Q9NP92": 10869,
    "P82664": 11070,
    "Q9BYC9": 11219,
    "Q9NRX2": 11394,
    "Q7Z7H8": 11655,
    "Q9BQ48": 11747,
    "Q7Z7F7": 11875,
    "P82650": 12235,
    "P82663": 12408,
    "Q9Y3D3": 12545,
    "Q9BYN8": 12750,
    "Q9BYD6": 13075,
    "Q92552": 13489,
    "P82912": 13683,
    "Q9Y2R5": 13813,
    "Q8N5N7": 13971,
    "Q6P1L8": 14116,
    "Q9BYD2": 14383,
    "Q96GC5": 14595,
    "P49406": 14887,
    "Q9H2W6": 15166,
    "Q9NX20": 15417,
    "P51398": 15815,
    "Q5T653": 16120,
    "Q9P0J6": 16223,
    "Q9HD33": 16473,
    "Q13405": 16639,
    "O75394": 16704,
    "Q4U2R6": 16832,
    "Q86TS9": 16955,
    "P82673": 17278,
    "Q96EL2": 17445
}

df = pd.read_csv("../bio_projects_files/coroletion_mutaion_result.csv")
i_protein = df['Residue_i']
j_protein = df['Residue_j']
for i, j in zip(i_protein, j_protein):
    i_number = int(i.split()[1])
    protein_i = get_protein_by_place(i_number, human_refrence)

    j_number = int(j.split()[1])
    protein_j = get_protein_by_place(j_number, human_refrence)
    #TODO: check how to do it - treshold to high
    if protein_j == protein_i:
        continue
    counting_matrix[protein_i][protein_j] += 1

# plt.imshow(counting_matrix, cmap='hot', interpolation='nearest')
# plt.show()

fig = px.imshow(counting_matrix, color_continuous_scale='RdBu_r', origin='lower')
# fig = px.imshow(counting_matrix, color_map='jet')
fig.show()