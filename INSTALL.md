# Installing PAM Scanning (easy, no coding)

This guide is for lab members and other outside researchers who just want to
**run** PAM Scanning. You do not
need to know any programming. There are only two steps, and the second is a
one-time setup that happens automatically.

PAM Scanning needs a free tool called **NCBI BLAST+**. To keep things simple, we
install it together with the app using **Miniforge** (a small, free program that
manages scientific software). You install Miniforge once; after that everything
else is automatic.

---

## Step 1 — Install Miniforge (one time, ~5 minutes)

1. Go to <https://conda-forge.org/download/>.
2. Download the installer for your computer:
   - **Mac**: choose the macOS installer (Apple Silicon for M1/M2/M3/M4, or Intel
     for older Macs — the page detects this for you).
   - **Windows**: choose the Windows installer.
3. Run the installer and accept the default options. (On Mac it opens in a
   Terminal window; just follow the prompts and let it finish.)
4. **Restart** your Terminal (Mac) or any open Command Prompt (Windows) so the
   install takes effect. On Mac you can also just restart the computer.

> You only ever do Step 1 once per computer.

## Step 2 — Get PAM Scanning and double-click to run

1. Get the code. pam_scanning is a **public** repository, so no account or password
   is needed — pick whichever way you prefer:
   - **Download ZIP (simplest):** on the [GitHub page](https://github.com/isomlab/pam_scanning),
     click the green **Code** button → **Download ZIP**, then unzip it somewhere easy
     to find, like your Desktop or Documents.
   - **GitHub Desktop (best if you'll update often):** [desktop.github.com](https://desktop.github.com)
     → **File ▸ Clone repository… ▸ URL** → paste
     `https://github.com/isomlab/pam_scanning` → **Clone**. Update later with
     **Fetch/Pull origin**.
   - **Terminal:** `git clone https://github.com/isomlab/pam_scanning.git`
     (update later with `git pull`).
2. Open the unzipped folder, then open the **`launchers`** folder inside it.
3. Double-click the launcher for your computer:
   - **Mac**: **`PAM Scanning.command`**
   - **Windows**: **`PAM Scanning.bat`**

The **first time** you double-click it, a window opens and sets things up
(it downloads BLAST+ and the app — this can take a few minutes). When that
finishes, the PAM Scanning window appears.

**Every time after that**, double-clicking the launcher just opens the app.

---

## Notes for Mac users

The first time you open `PAM Scanning.command`, macOS may warn that it is from
an unidentified developer. To allow it:

- **Right-click** (or Control-click) the file → **Open** → **Open** in the dialog.

You only need to do this the first time. (This is normal macOS security for files
downloaded from the internet.)

## If something goes wrong

- **"Could not find conda"** — Miniforge isn't installed yet, or the Terminal /
  Command Prompt wasn't restarted after installing it. Do Step 1 again and
  restart, then try the launcher.
- **Setup didn't finish** — make sure you are connected to the internet (the
  first run downloads BLAST+), then double-click the launcher again.
- Still stuck? Send the lab the message shown in the launcher window.

---

## For people comfortable with a terminal

You don't need the launchers — from the project folder:

```bash
conda env create -f environment.yml   # one-time setup (installs BLAST+ too)
conda activate pam_scanning
pam-scan-gui                           # graphical app
pam-scan --help                        # command-line version
```

See [`docs/usage.md`](docs/usage.md) for all options and
[`docs/blast_setup.md`](docs/blast_setup.md) for building a BLAST database.

### "No module named 'pam_scanning.gui'"

If a command fails with:

```
ModuleNotFoundError: No module named 'pam_scanning.gui'
```

another copy of `pam_scanning` is shadowing the installed one. Directories on
`PYTHONPATH` take precedence over the environment's `site-packages`, and so does the
current working directory, so an older copy of this project in either place wins.

The double-click launchers already guard against this. If you are running the commands
yourself, check what actually got imported:

```bash
python -c "import pam_scanning; print(pam_scanning.__file__)"
```

If that path is not inside your environment's `site-packages`, run with `PYTHONPATH`
cleared:

```bash
env -u PYTHONPATH pam-scan-gui        # macOS / Linux
```

```bat
set PYTHONPATH=  &  pam-scan-gui      REM Windows
```

To fix it permanently, rename the older copy or drop its parent directory from
`PYTHONPATH`.
