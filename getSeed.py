import random

master_seed = 35       # whatever --seed you used
n_cells     = 88     # whatever --n-cells you used
cell_number = 4        # the cell you want

master_rng  = random.Random(master_seed)
cell_seeds  = [master_rng.randint(0, 2**31 - 1) for _ in range(n_cells)]

print(f"Seed for cell {cell_number}: {cell_seeds[cell_number - 1]}")