import pandas as pd
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

df = pd.read_excel(r'Peptide complex dataset.xlsx',header=1)

print(df)

df = df[df['Interface score']>= -1000]
score = df['Interface score']

mu, std = norm.fit(score)
plt.hist(score, bins=25, density=True, alpha=0.6, color='g')

mu = -49.48
std =25.52
interface_score_range = [mu-std,mu+std]