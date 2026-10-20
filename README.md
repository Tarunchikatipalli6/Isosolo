# isosolo — Getting Started

A complete guide for someone who has just downloaded this folder and never run it before.

---

## Before you begin

**What you need:**

- A computer running macOS, Linux, or Windows with WSL
- An internet connection (the tool fetches sequences from Ensembl)
- About 15 minutes for first-time setup

**What you don't need:**

- Any prior Python knowledge
- A high-performance cluster — a laptop is fine
- Your own sequence files

**What this tool does:** You give it a gene name and one transcript ID. It gives you a PCR primer pair that amplifies that transcript and not the other versions of the same gene.

---

## Part 1 — Setup (once)

### 1.1 Open a terminal

- **macOS:** press `Cmd + Space`, type "Terminal", press Enter
- **Linux:** `Ctrl + Alt + T`
- **Windows:** open "Ubuntu" or "WSL" from the Start menu

A window with a text prompt appears. Every command below is typed there, one line at a time, pressing Enter after each.

### 1.2 Go to the downloaded folder

If you downloaded a ZIP from GitHub, unzip it first by double-clicking. Then:

```bash
cd ~/Downloads/isosolo
```

Adjust the path if you put it elsewhere. Check you're in the right place:

```bash
ls
```

You should see `pyproject.toml`, `README.md`, and a folder called `src`. If you see something else, you're in the wrong directory, or the real folder is nested one level deeper — try `cd isosolo` again.

### 1.3 Check you have Python

```bash
python3 --version
```

You need **3.9 or newer**. If you get "command not found," install Python from [python.org](https://www.python.org/downloads/) first.

### 1.4 Install

```bash
pip install -e .
```

This reads `pyproject.toml` and pulls in the two things isosolo needs: `primer3-py` for melting-temperature calculations and `requests` for fetching sequences.

Takes a minute or two. Lots of text scrolls past — normal.

**If `pip` isn't found:** try `python3 -m pip install -e .`

**If you get a permissions error:** use a virtual environment instead (recommended anyway):

```bash
python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -e .
```

With a venv you must run `source venv/bin/activate` each time you open a new terminal.

### 1.5 Confirm it worked

```bash
isosolo --version
```
