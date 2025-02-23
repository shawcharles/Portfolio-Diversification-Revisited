from scipy.io import loadmat
import numpy as np
import os

# Get the current directory and construct the full path
current_dir = os.path.dirname(os.path.abspath(__file__))
matlab_file = os.path.join(current_dir, 'matlab_results.mat')

print(f"Loading MATLAB file from: {matlab_file}")
data = loadmat(matlab_file)

# Print information about DataState1
print("\nDataState1 info:")
if 'DataState1' in data:
    state1 = data['DataState1']
    print(f"Shape: {state1.shape}")
    print(f"Unique values: {np.unique(state1)}")
else:
    print("DataState1 not found in file")

# Print information about DataState2
print("\nDataState2 info:")
if 'DataState2' in data:
    state2 = data['DataState2']
    print(f"Shape: {state2.shape}")
    print(f"Unique values: {np.unique(state2)}")
else:
    print("DataState2 not found in file")

# Print information about Data
print("\nData info:")
if 'Data' in data:
    returns = data['Data']
    print(f"Shape: {returns.shape}")
    print(f"Min: {np.min(returns)}")
    print(f"Max: {np.max(returns)}")
    print(f"Mean: {np.mean(returns)}")
else:
    print("Data not found in file")

# Print information about P_Rt_1 and P_Rt_2 (regime probabilities)
print("\nRegime Probabilities info:")
if 'P_Rt_1' in data and 'P_Rt_2' in data:
    p_rt_1 = data['P_Rt_1']
    p_rt_2 = data['P_Rt_2']
    print("P_Rt_1:")
    print(f"Shape: {p_rt_1.shape}")
    print(f"Min: {np.min(p_rt_1)}")
    print(f"Max: {np.max(p_rt_1)}")
    print("\nP_Rt_2:")
    print(f"Shape: {p_rt_2.shape}")
    print(f"Min: {np.min(p_rt_2)}")
    print(f"Max: {np.max(p_rt_2)}")
else:
    print("Regime probabilities not found in file")
