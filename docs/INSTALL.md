# Install from scratch (first time — no conda yet)

This is the **one-time setup of your computer** so it can run PAM Scanning. It assumes
you have installed *nothing* yet. You don't need to know how to code — do these
two things once, in order (about 15 minutes). When you're done, go to
**[getting_started.md](getting_started.md)** to download and run the app.

> **Mac vs Windows:** steps are the same. Where Mac says **Terminal**, Windows users
> use the **Miniforge Prompt** app (installed below). Extra Windows notes are called
> out as they come up.

> **No GitHub account needed.** This is a public repository — anyone can
> download it without signing up for anything.

---

## 1. Install Miniforge (free)

Miniforge gives you a private copy of Python plus everything the app needs, without
disturbing anything else on your computer.

**Mac:**
1. Go to the Miniforge download page:
   **[conda-forge.org/download](https://conda-forge.org/download/)**.
2. Choose the **macOS** installer that matches your Mac:
   - **Apple Silicon** (`arm64`) — Macs from ~2020 on. *(Not sure? Apple menu  →
     **About This Mac**. If it says "Chip: Apple M1/M2/M3/M4", you're Apple Silicon.)*
   - **Intel** (`x86_64`) — older Macs ("Intel Core i5/i7").
   - Use the **`.pkg`** (graphical) installer — it's the click-through one.
3. Double-click the downloaded `.pkg` and click through, **accepting the defaults**.
4. **Quit Terminal if it's open, then open a fresh Terminal** (press ⌘-Space, type
   `Terminal`, press Enter). If you see **`(base)`** at the start of the line, conda
   is installed correctly.

**Windows:**
1. Download the **Windows** Miniforge installer from the same page and run it,
   accepting the defaults.
2. From the Start menu, open **Miniforge Prompt** — that's the window you'll type in
   for the next guide. (You should see `(base)` at the start of the line.)

---

## 2. Install GitHub Desktop (optional, free)

This lets you download and later update the code with buttons instead of typing.

1. Go to **[desktop.github.com](https://desktop.github.com)**, download, and install.
2. Open it. You can **skip the sign-in** — this repo is public, so you don't need an account to clone it.

*(You can skip this and download a ZIP instead — the next guide shows both.)*

---

## You're set up

Your computer now has:

- ✅ **conda** (Python + the tools the app needs), and
- ✅ **GitHub Desktop** (optional, for easy downloads and updates).

You only do the above **once**. Next, follow
**[getting_started.md](getting_started.md)** to download the code and launch
PAM Scanning — that part is a double-click.

---

### Trouble?

- **After installing Miniforge you don't see `(base)`** — fully quit Terminal and open
  a brand-new window; the installer only takes effect in a fresh one. On Mac, if it
  still doesn't show, run `source ~/miniforge3/bin/activate` once.
- **Which Mac installer?** Apple menu  → About This Mac → the "Chip"/"Processor"
  line tells you Apple Silicon (M-series) vs Intel.
- **Already have Anaconda or Miniconda?** That works too — you don't need Miniforge as
  well. Any `conda` on your PATH is fine.

