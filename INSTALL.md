# PORTEK Installation Guide

## Requirements

PORTEK requires **Python 3.12 or higher**.
Make sure you are in PORTEK directory when creating the virtual environment and installing the package.

### Checking your Python version
```bash
python --version
```

If you don't have Python 3.12, consider using [pyenv](https://github.com/pyenv/pyenv) to install and manage Python versions:

```bash
# Install Python 3.12 with pyenv (optional)
pyenv install 3.12
pyenv local 3.12
```

## Installation in a Virtual Environment

### 1. Create a virtual environment
```bash
python3.12 -m venv env
# Or if pyenv is set to 3.12:
python -m venv env
```

### 2. Activate the virtual environment
```bash
source env/bin/activate  # On Linux/Mac
# OR
env\Scripts\activate  # On Windows
```

### 3. Install PORTEK in development mode
```bash
pip install -e .
```

This will install PORTEK and all its dependencies with pinned versions.

### 4. (Optional) Install testing dependencies
```bash
pip install -e ".[test]"
```

## Usage

After installation, you can run PORTEK using the `portek` command instead of `python PORTEKrun.py`:

```bash
# Create a new project
portek new $project_directory

# Find optimal k value
portek find_k $project_directory --min_k 5 --max_k 31

# Find enriched k-mers
portek find_enriched $project_directory -k 15

# Map k-mers to reference
portek map $project_directory -k 15 -d 2

# Construct phylogenetic tree
portek tree $project_directory -k 15
```

## Deactivating the Environment

When you're done working:
```bash
deactivate
```

## Uninstalling

To uninstall PORTEK:
```bash
pip uninstall portek
```
