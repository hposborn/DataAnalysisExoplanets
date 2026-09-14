---
title: Installing the necessary code & packages
---

# Git

Git is version control system. It allows you to track changes in your files.

_NB: *Git* is different from *GitHub* which is simply a respository to store version controlled code online_ [+]

You can interact with it **in the terminal** (i.e. unix terminal, or in e.g. `zed`/`vsc`)

-v-

## Downloading & syncing this courses package

In a terminal, navigate to where you want to download the package and run:

`git clone http://github.com/hposborn/DataAnalysisForExoplanets`

_(You may need to install git on your system first with e.g. `apt-get git`, `brew pour git`, etc)_
<!-- .element: class="fragment" -->

-v-


## Making (and saving) your local changes

You will be modifying and playing with the python files provided here. In order to save these files, it is best to frequently **commit** your changes. However, you don't want to overwrite what everyone else is looking at, so you should create a **branch**.

### Creating a branch:
`git checkout -b my-changes`
This swaps from the `main` branch to your own personal branch

### Commiting changes
Now you can commit the code (and add comments):
`git commit -am "Finished Gaia/Stars project"`

-v-

### Updating the main branch
Every week, to ensure changes I may appear in your project folder, you need to `pull` the `main` branch. This requires:
1) Ensure your personal branch is commited using `git commit`
<!-- .element: class="fragment" -->
2) Swapping back onto the main branch using `git checkout main`
<!-- .element: class="fragment" -->
3) Pulling the new main code using `git pull`
<!-- .element: class="fragment" -->
4) Then go back into your personal branch (`git checkout my-changes`)
<!-- .element: class="fragment" -->
5) Merge the local branch with the main using `git merge main` 
<!-- .element: class="fragment" -->
6) Resolve conflicts (e.g. using VSC _Source Control_ tab) 
<!-- .element: class="fragment" -->

---

# Zed

This is a much quicker development environment (thanks to its Rust back-end) which does not cook your laptop (unlike VS Code).

Download at https://zed.dev/download; or on Mac, you can use `brew install zed`.

-v-

We will use REPL (read-evalutate-print-loop; see https://zed.dev/docs/repl)

Allows .py files to be "inspected" like notebooks using the `# %%` command.

The command `# %%` is used to separate cells.
<!-- .element: class="fragment" -->

Outputs (e.g. plots) can be displayed just like a notebooks with VS Code.
<!-- .element: class="fragment" -->

To run each cell in Zed press Shift + Control + Enter.
<!-- .element: class="fragment" -->
-v-


### Python kernels

You will need to ensure you are running the correct python environment. On Zed/VSC, this can be modified using the top-right button. 

But first let's set up the environment with `uv`...

---

# UV python environment

## Installation:

### On Linux/Mac:
`curl -LsSf https://astral.sh/uv/install.sh | sh`

### On Windows:
`powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`

-v-

# Set-up

### All you need to run:

`uv sync`

Make sure to run this every new project (as I may add additional packages here)

This environment can then be activated using `source .venv/bin/activate`

### Installing packages

You can then use `uv add` to add python modules.

-v-


## To create a new python environment

`uv venv --python 3.14` or more simply `uv run` (e.g. when .venv file exists)

This can then be activated using `source .venv/bin/activate`

You can also directly run python scripts without launching the environment via `uv run main.py`

-v-

## Running REPL in Zed

You may need to run `python -m ipykernel install --user` to ensure the repl functionality works.

---

# Python installation issues?

### For Mac
- You will likely need to download XCode (10+Gb)

### For linux
- Should be native...

### For windows
- No idea (sorry!)

---

# Viewing the slides

To "build" the slides:
1) Make sure to have the most up-to-date version using `git pull`
2) Build the slides with `mkslides build all_slides/`
3) Set-up the interactive slide back-end `mkslides serve all_slides/`
4) Open `http://localhost:8000` in a browser
5) View the slides in HTML, or download as PDF.
