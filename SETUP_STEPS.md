# AmberFlow_development Repository Setup Steps

This document outlines all the steps taken to create the `AmberFlow_development` repository and sync files from the original `AmberFlow` repository.

## Steps Performed

### 1. Create New Folder
```bash
mkdir -p AmberFlow_development
```
Created a new directory called `AmberFlow_development` in the home directory (`/home/hn533621/`).

### 2. Copy All Files from AmberFlow
```bash
cp -r AmberFlow/* AmberFlow_development/
cp -r AmberFlow/.[!.]* AmberFlow_development/
```
Copied all files and hidden files (except `.` and `..`) from the `AmberFlow` directory to the new `AmberFlow_development` directory.

**Files copied:**
- `add_caps.py`
- `css/` directory
- `Dockerfile`
- `.gitattributes`
- `html/` directory
- `js/` directory
- `obsolete/` directory
- `output/` directory
- `python/` directory
- `README.md`
- `start_web_server.py`

### 3. Remove Old Git Repository
```bash
rm -rf .git
```
Removed the existing `.git` directory from the copied files to start fresh with a new git repository.

### 4. Create .gitignore File
Created a `.gitignore` file to exclude unnecessary files from version control:

**Contents:**
- Python cache files (`__pycache__/`, `*.pyc`, etc.)
- Output files (`output/`, `*.pdb`, `*.cif`, `*.log`)
- IDE files (`.vscode/`, `.idea/`, etc.)
- OS files (`.DS_Store`, `Thumbs.db`)
- Environment files (`.env`, `venv/`, etc.)

### 5. Initialize Git Repository
```bash
git init
```
Initialized a new empty git repository in the `AmberFlow_development` directory.

### 6. Rename Branch to main
```bash
git branch -m main
```
Renamed the default branch from `master` to `main` to follow modern git conventions.

### 7. Stage All Files
```bash
git add .
```
Staged all files in the repository for the initial commit.

### 8. Create Initial Commit
```bash
git commit -m "Initial commit: Copy from AmberFlow"
```
Created the initial commit with all copied files.

**Commit details:**
- 11 files changed
- 6,567 insertions
- Files included: `.gitattributes`, `.gitignore`, `Dockerfile`, `README.md`, `add_caps.py`, CSS, HTML, JS files, Python application files, and `start_web_server.py`

### 9. Create Private GitHub Repository
```bash
gh repo create AmberFlow_development --private --source=. --remote=origin --push
```
Created a private GitHub repository named `AmberFlow_development` and pushed all files.

**Repository details:**
- Repository URL: https://github.com/nagarh/AmberFlow_development
- Visibility: Private
- Remote name: `origin`
- Branch: `main`

### 10. Verify Setup
```bash
git remote -v
git status
```
Verified that:
- Remote `origin` is correctly configured to point to the GitHub repository
- All files are committed and pushed
- Working tree is clean

## Summary

The `AmberFlow_development` repository has been successfully created as a private GitHub repository with all files from the original `AmberFlow` directory. The repository is now ready for development work and version control.

## Future Git Operations

To make changes and push them to GitHub:

1. **Stage changes:**
   ```bash
   git add .
   ```

2. **Commit changes:**
   ```bash
   git commit -m "Your commit message"
   ```

3. **Push to GitHub:**
   ```bash
   git push
   ```

## Notes

- The original `AmberFlow` repository remains unchanged and still points to the Hugging Face Space
- The new `AmberFlow_development` repository is completely independent
- A `.gitignore` file has been added to exclude unnecessary files from version control

