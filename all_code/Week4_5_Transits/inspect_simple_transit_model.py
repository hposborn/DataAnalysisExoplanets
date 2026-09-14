# %% [markdown]

# # Let's build a very simple transit model!

# First load some modules...

# %% Importing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %%

t=np.linspace(0,1,100)
y=np.sin(t)+np.random.random(len(t))
print("this output")
plt.plot(t,y)
plt.show()
# %%
