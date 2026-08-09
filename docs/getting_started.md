# Getting started with PAM Scanning

PAM Scanning designs CRISPR/Cas9 guide RNAs and chimera-insertion primers across an
ORF, with synonymous PAM silencing and BLAST+ off-target screening.

Everything the lab tools have in common — downloading, launching, updating, and what
to do when something goes wrong — is on one shared page: **[Getting started with a lab
tool](https://dangerisom.github.io/Isom-Lab/getting-started/)**. This guide covers
only what is specific to PAM Scanning.

---

## Before you start

Your computer needs **conda** (Miniforge, Miniconda, or Anaconda).

> **First time on this computer?** Do the one-time **[Setting up your
> computer](https://dangerisom.github.io/Isom-Lab/setup/)** first, then come back here.
> PAM Scanning's own install notes are in **[INSTALL.md](INSTALL.md)**.

---

## Get it and launch it

**1. Download it** from
   **[github.com/isomlab/pam_scanning](https://github.com/isomlab/pam_scanning)**. It
   is **public**, so no account or password is needed: **Download ZIP**, **GitHub
   Desktop**, or `git clone`. Step by step: **[Get the
   code](https://dangerisom.github.io/Isom-Lab/getting-started/#public-tools)**.

**2. Open the `launchers` folder inside it and double-click:**

- **Mac:** `PAM Scanning.command`
- **Windows:** `PAM Scanning.bat`

The **first** launch takes a few minutes while it builds a private, isolated conda
environment named `pam_scanning` containing Python and everything the app needs. Every
launch after that opens straight away. You don't need to type anything.

New to this? **[Launch
it](https://dangerisom.github.io/Isom-Lab/getting-started/#launch-it)** walks through
what you will see, and **[the first-time
hiccups](https://dangerisom.github.io/Isom-Lab/getting-started/#the-first-time-hiccups)**
covers macOS blocking the file, Windows SmartScreen, and a double-click that does
nothing.

---

## Use it

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

Refresh the folder the way you got it — see **[Updating
later](https://dangerisom.github.io/Isom-Lab/getting-started/#updating-later)**. If a
release changes what the tool depends on, delete its environment and let the launcher
rebuild it on the next double-click:

```bash
conda env remove -n pam_scanning
```

---

## If something goes wrong

Most problems are one of a handful of things — start with **[Try these
first](https://dangerisom.github.io/Isom-Lab/getting-started/#try-these-first)**, then
the rest of **[If something goes
wrong](https://dangerisom.github.io/Isom-Lab/getting-started/#if-something-goes-wrong)**.

Stuck? Send Dan (<a href="&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#100;&#105;&#115;&#111;&#109;&#64;&#109;&#105;&#97;&#109;&#105;&#46;&#101;&#100;&#117;">&#100;&#105;&#115;&#111;&#109;<span>&#64;</span>&#109;&#105;&#97;&#109;&#105;<span>&#46;</span>&#101;&#100;&#117;</a>) the exact command you ran and the message you got — the shared page lists **[what to include](https://dangerisom.github.io/Isom-Lab/getting-started/#still-stuck)**.
