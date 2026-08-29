"""Tests for the GUI's last-directory persistence (no display required).

Importing :mod:`pam_scanning.gui` only needs the tkinter module to import, not a
running display, so these run headless. Only the pure state helpers are exercised.
"""

import os

from pam_scanning import gui


def test_last_dir_round_trips(tmp_path, monkeypatch):
    state = tmp_path / "nested" / "gui_state.json"
    monkeypatch.setattr(gui, "_STATE_PATH", state)
    gui._save_last_dir(str(tmp_path))
    assert gui._load_last_dir() == str(tmp_path)
    assert state.is_file()  # parent directory was created


def test_missing_state_falls_back_to_cwd(tmp_path, monkeypatch):
    monkeypatch.setattr(gui, "_STATE_PATH", tmp_path / "absent.json")
    assert gui._load_last_dir() == os.getcwd()


def test_stale_directory_is_ignored(tmp_path, monkeypatch):
    monkeypatch.setattr(gui, "_STATE_PATH", tmp_path / "gui_state.json")
    gui._save_last_dir(str(tmp_path / "was_deleted"))  # never existed
    assert gui._load_last_dir() == os.getcwd()


def test_corrupt_state_is_tolerated(tmp_path, monkeypatch):
    state = tmp_path / "gui_state.json"
    state.write_text("{ not valid json")
    monkeypatch.setattr(gui, "_STATE_PATH", state)
    assert gui._load_last_dir() == os.getcwd()  # no exception


# --- scroll direction ----------------------------------------------------
# The sign here was got wrong twice by hard-coding a reversal. These lock it to
# Tk's own convention: ::tk::MouseWheel divides by a NEGATIVE factor, so a
# positive delta scrolls the view UP, and X11 button 4 is wheel-up.

def test_positive_delta_scrolls_the_view_up():
    from pam_scanning.gui import wheel_step
    assert wheel_step(None, 120) == -1
    assert wheel_step(None, 3) == -1


def test_negative_delta_scrolls_the_view_down():
    from pam_scanning.gui import wheel_step
    assert wheel_step(None, -120) == 1
    assert wheel_step(None, -3) == 1


def test_x11_buttons_follow_the_same_convention():
    from pam_scanning.gui import wheel_step
    assert wheel_step(4, 0) == -1     # button 4 is wheel up
    assert wheel_step(5, 0) == 1      # button 5 is wheel down


def test_zero_delta_does_not_scroll():
    from pam_scanning.gui import wheel_step
    assert wheel_step(None, 0) == 0


def test_direction_matches_tk_mousewheel_for_a_range_of_deltas():
    """tk::MouseWheel is `yview scroll [expr {$amount/$factor}]` with factor<0."""
    from pam_scanning.gui import wheel_step
    for delta in (-240, -120, -40, -1, 1, 40, 120, 240):
        tk_sign = -1 if (delta / -40.0) < 0 else 1
        assert wheel_step(None, delta) == tk_sign, delta
