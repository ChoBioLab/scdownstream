#!/usr/bin/env python3

import platform
import os
import pickle

import anndata as ad
import pandas as pd
import numpy as np
import yaml

adata = ad.read_h5ad("${adata}")
prefix = "${prefix}"
obs_paths = sorted("${obs_paths}".split())
var_paths = sorted("${var_paths}".split())
obsm_paths = sorted("${obsm_paths}".split())
obsp_paths = sorted("${obsp_paths}".split())
uns_paths = sorted("${uns_paths}".split())
layers_paths = sorted("${layers_paths}".split())

def simple_name(path):
    basename = os.path.basename(path)
    return basename[:basename.rfind(".")]

def load_pickle_or_csv(path):
    if path.endswith(".pkl"):
        data = pd.read_pickle(path)
    elif path.endswith(".csv"):
        data = pd.read_csv(path, index_col=0)
    else:
        raise ValueError(f"Unsupported file extension: {path}")

    # Ensure we return a DataFrame
    if isinstance(data, pd.Series):
        data = data.to_frame()

    return data

def fix_dtypes(data):
    """Fix problematic dtypes for HDF5 compatibility"""
    # Handle Series (single column)
    if isinstance(data, pd.Series):
        if data.dtype == 'object':
            # Handle boolean columns that got mixed types
            if 'highly_variable' in data.name or 'variable' in data.name:
                return data.astype(bool)
            # Handle string columns
            elif data.apply(lambda x: isinstance(x, str) if pd.notna(x) else True).all():
                return data.astype(str)
            # Try numeric conversion
            else:
                try:
                    numeric_data = pd.to_numeric(data, errors='coerce')
                    if numeric_data.isna().all():
                        return data.astype(str)
                    return numeric_data
                except:
                    return data.astype(str)
        return data

    # Handle DataFrame (multiple columns)
    elif isinstance(data, pd.DataFrame):
        for col in data.columns:
            if data[col].dtype == 'object':
                # Handle boolean columns that got mixed types
                if 'highly_variable' in col or 'variable' in col:
                    data[col] = data[col].astype(bool)
                # Handle string columns
                elif data[col].apply(lambda x: isinstance(x, str) if pd.notna(x) else True).all():
                    data[col] = data[col].astype(str)
                # Try numeric conversion
                else:
                    try:
                        data[col] = pd.to_numeric(data[col], errors='coerce')
                        if data[col].isna().all():
                            data[col] = data[col].astype(str)
                    except:
                        data[col] = data[col].astype(str)
        return data

    # If neither Series nor DataFrame, return as-is
    return data

for path in obs_paths:
    if path:  # Check if path is not empty
        df = load_pickle_or_csv(path).reindex(adata.obs_names)
        df = fix_dtypes(df)
        adata.obs = pd.concat([adata.obs, df], axis=1)

for path in var_paths:
    if path:  # Check if path is not empty
        df = load_pickle_or_csv(path).reindex(adata.var_names)
        df = fix_dtypes(df)
        adata.var = pd.concat([adata.var, df], axis=1)

# Final dtype cleanup
adata.var = fix_dtypes(adata.var)
adata.obs = fix_dtypes(adata.obs)

for path in obsm_paths:
    if path:  # Check if path is not empty
        df = pd.read_pickle(path).reindex(adata.obs_names)
        adata.obsm[simple_name(path)] = np.float32(df.to_numpy())

for path in obsp_paths:
    if path:  # Check if path is not empty
        adata.obsp[simple_name(path)] = np.load(path, allow_pickle=True).item()

for path in uns_paths:
    if path:  # Check if path is not empty
        adata.uns[simple_name(path)] = pickle.load(open(path, "rb"))

for path in layers_paths:
    if path:  # Check if path is not empty
        adata.layers[simple_name(path)] = np.float32(np.load(path))

adata.write_h5ad(f"{prefix}.h5ad")
adata.obs.to_csv(f"{prefix}_metadata.csv")

# Versions
versions = {
    "NFCORE_SCDOWNSTREAM:SCDOWNSTREAM:FINALIZE:ADATA_EXTEND": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
