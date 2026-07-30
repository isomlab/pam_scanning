# Getting started with PAM Scanning

This guide gets PAM Scanning running on your computer. You do **not** need to know how
to code — follow the steps and copy-paste where asked. It takes about 10 minutes,
once. After that, launching it is a double-click.

> On Windows everything below works the same — use the **Miniforge Prompt** app
> wherever this says "Terminal", and see the Windows note in each step.

---

## Before you start

This guide assumes your computer already has **conda** (Miniforge, Miniconda, or
Anaconda).

> **First time? Never installed conda?** Do the one-time
> **[install-from-scratch guide → INSTALL.md](INSTALL.md)** first, then come back here.

---

## Step 1 — Get the code

PAM Scanning is a **public** repository, so no account or password is needed. Pick
whichever way you prefer.

**Option A — Download ZIP (fastest, nothing to install).**
1. Open **[github.com/isomlab/pam_scanning](https://github.com/isomlab/pam_scanning)**.
2. Click the green **`Code ▾`** button → **Download ZIP**.
3. Double-click the downloaded file to unzip it. You now have a folder called
   **`pam_scanning-main`** — move it somewhere easy, like your **Documents**.

**Option B — GitHub Desktop (best if you'll update often).**
1. Open GitHub Desktop → **File ▸ Clone repository… ▸ URL**.
2. Paste `https://github.com/isomlab/pam_scanning` → pick a folder (e.g. Documents) → **Clone**. This makes a folder
   called **`pam_scanning`**.

**Option C — `git clone` in Terminal.** Because the repo is public this just works,
with no password:

```bash
cd ~/Documents
git clone https://github.com/isomlab/pam_scanning.git
```

Either way you now have a pam_scanning folder on your computer.

---

## Step 2 — Launch it

**Open the `launchers` folder inside that folder and double-click:**

- **Mac:** `PAM Scanning.command`
- **Windows:** `PAM Scanning.bat`

That's it. The **first** launch takes a few minutes: it builds a private, isolated
conda environment (named `pam_scanning`) containing Python, the app, and the external **NCBI BLAST+** toolkit it depends on, then opens the app.
**Every launch after that opens straight away.**

You don't need to type anything, and it won't touch any other Python on your computer.

> **Mac note:** the first time, macOS may say the file is from an unidentified
> developer. Right-click (or Control-click) the file → **Open** → **Open**. You only
> do this once.

> **Mac, if double-clicking does nothing:** the file may have lost its executable
> flag when unzipped. In Terminal, run
> `chmod +x "<your folder>/launchers/PAM Scanning.command"` once, then double-click again.

---

## Step 3 — Use it

1. **Choose your ORF** and its flanking genomic sequence.
2. **Pick the codons** to scan — every codon by default, or a chosen subset.
3. **Run.** The app finds PAM sites, silences them with synonymous substitutions,
   screens guides for off-targets by BLAST, picks the best guide per insertion
   codon, and writes primer orders, plate layouts, and QC reports.

Off-target screening needs a local BLAST database — see
[`docs/blast_setup.md`](blast_setup.md). Full pipeline description:
[`docs/pipeline.md`](pipeline.md); day-to-day options: [`docs/usage.md`](usage.md).

If you'd rather work in a Terminal:

```bash
conda activate pam_scanning
pam-scan --help              # all options
pam-scan-fetch-cds --help    # fetch a CDS from UniProt
```

---

## Updating later

- **GitHub Desktop (Option B):** open it and click **Fetch / Pull origin**.
- **`git clone` (Option C):** `cd ~/Documents/pam_scanning && git pull`.
- **Downloaded the ZIP (Option A):** download a fresh ZIP and replace the old
  folder's contents (keep the same folder name and location).

The environment installs the code in "editable" mode, so an update takes effect the
next time you launch — no reinstall. If a release changes the dependencies, delete the
environment and let the launcher rebuild it:

```bash
conda env remove -n pam_scanning
```

---

## If something goes wrong

- **"conda: command not found"** — close and reopen Terminal after installing conda
  (the installer needs a fresh window). On Mac, if it still isn't found, run
  `source ~/miniforge3/bin/activate` once.
- **The launcher says it can't find conda** — same cause. Install conda from the
  [install-from-scratch guide](INSTALL.md), then double-click the launcher again.
- **"pam-scan-gui: command not found"** — you probably forgot
  `conda activate pam_scanning` first. Run it, then try again.
- **The window doesn't appear** — make sure you used the **conda** install above; its
  Python includes the graphics toolkit the app needs. A plain system-Python
  `pip install` can be missing it.
- **Setup failed partway through** — remove the half-built environment with
  `conda env remove -n pam_scanning` and double-click the launcher again.

Stuck? Send Dan the exact command you ran and the message you got.

---

## Alternative: plain `pip` (if you don't use conda)

Conda is recommended because it guarantees the GUI toolkit is present. If you'd rather
use `pip`, from inside the pam_scanning folder:

```bash
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e .
pam-scan-gui
```

This only works if your Python includes **tkinter**: macOS's built-in `python3` does;
Homebrew Python needs `brew install python-tk`; conda always does. Note that `pip`
does **not** install BLAST+ — see [`blast_setup.md`](blast_setup.md).
