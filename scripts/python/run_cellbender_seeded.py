#!/usr/bin/env python
"""
Reproducible CellBender wrapper with seed control
Runs CellBender remove-background with all random seeds set for reproducibility
"""
import os
import random
import sys
import numpy as np
import torch
from cellbender.base_cli import main # this is how we run cellbender from inside python

# Set seed for reproducibility
seed = 281330800

# --- Set All Random Seeds ---
print(f"Setting seeds to {seed} for reproducibility...")
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)

if torch.cuda.is_available():
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # for reproducibility on GPU
    
    # Our seeds may be overwritten by cellbender's constant seeds but
    # these cudnn settings are not and are essential for reproducibility. 
    # Without these, cudnn's non-deterministic algorithms can cause
    # variation between runs even with fixed random seeds.
    torch.backends.cudnn.deterministic = True 
    torch.backends.cudnn.benchmark = False
    print("CUDA seeds set. CUDNN is in deterministic mode.")
else:
    print("CUDA not available. Running on CPU (will be reproducible).")

# Print the command-line arguments being used
print(f"Starting CellBender with arguments: {sys.argv[1:]}")

# Modify sys.argv to include the 'remove-background' command
# CellBender CLI expects: ['cellbender', 'remove-background', '--input', ...]
sys.argv = ['cellbender', 'remove-background'] + sys.argv[1:]

# Run cellbender's main() function
try:
    main()
    print("CellBender run completed successfully.")
except Exception as e:
    print(f"CellBender run failed: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)