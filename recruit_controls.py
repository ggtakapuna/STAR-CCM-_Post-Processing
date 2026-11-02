"""
================================================================================
  CFD Slicer Pipeline
  Author  : John Baton
  Project : UTFR
--------------------------------------------------------------------------------
  recruit_controls.py
  Shared Tk controls for Run and Compare recruit allocation.
  Displays input-based RAM estimates, resource limits and processing progress.
  Do not distribute without permission.
================================================================================
"""

import math
import os
import queue
import threading
import tkinter as tk
from tkinter import ttk, messagebox

import config as cfg
import recruits

BG = '#f5f5f5'
ACCENT = '#1976D2'
FONT = ('Segoe UI', 9)
FONT_B = ('Segoe UI', 9, 'bold')


class RecruitControls:
    _comparison = False

    def _init_recruits(self):
        self._running = False
        self._analyzing = False
        self._updating_allocation = False
        self._memory_profile = None
        self._events = queue.Queue()
        self._analysis_events = queue.Queue()

    def _build_recruit_controls(self, f, r):
        tk.Frame(f, bg=ACCENT, height=2).grid(row=r, column=0, columnspan=3,
                                            sticky='ew', padx=8, pady=(16, 0))
        r += 1
        tk.Label(f, text='Recruits & RAM Budget', font=FONT_B, bg=BG, fg=ACCENT).grid(
            row=r, column=0, sticky='w', padx=8, pady=(2, 4))
        r += 1
        self.v_parallel = tk.BooleanVar(value=True)
        self.parallel_check = tk.Checkbutton(f, text='Use multiple recruits (up to 8 planes at once)',
            variable=self.v_parallel, command=self._update_allocation, font=FONT, bg=BG)
        self.parallel_check.grid(row=r, column=0, columnspan=3, sticky='w', padx=8)
        r += 1
        resources = recruits.system_resources()
        ram_max = self._ram_max(resources)
        self.v_ram = tk.DoubleVar(value=min(ram_max, math.floor(
            recruits.default_budget_gib(resources)*4)/4))
        tk.Label(f, text='RAM budget (GiB)', font=FONT, bg=BG).grid(
            row=r, column=0, sticky='w', padx=8)
        self.ram_slider = tk.Scale(f, from_=0, to=ram_max, resolution=0.25,
            orient='horizontal', variable=self.v_ram, bg=BG, highlightthickness=0)
        self.ram_slider.grid(row=r, column=1, columnspan=2, sticky='ew', padx=8)
        r += 1
        self.allocation_lbl = tk.Label(f, text='', font=FONT_B, bg=BG, fg=ACCENT,
                                       justify='left', anchor='w', wraplength=650)
        self.allocation_lbl.grid(row=r, column=0, columnspan=3, sticky='ew', padx=8)
        r += 1
        tk.Label(f, text='Each recruit processes one plane and all selected scalars.\n'
                 'RAM is an estimate, not a hard limit. The count is recalculated after loading;\n'
                 'the parent process also needs memory to load the input files.',
                 font=FONT, bg=BG, fg='#555', justify='left', wraplength=650).grid(
                     row=r, column=0, columnspan=3, sticky='w', padx=8, pady=4)
        r += 1
        self.progress = ttk.Progressbar(f, mode='determinate')
        self.progress.grid(row=r, column=0, columnspan=3, sticky='ew', padx=8, pady=4)
        r += 1
        self.analyze_button = tk.Button(f, text='Analyze RAM from input files',
            command=self._analyze_ram, bg='#757575', fg='white', font=FONT,
            relief='flat', padx=10, pady=5, cursor='hand2')
        self.analyze_button.grid(row=r, column=0, columnspan=3, sticky='w', padx=8, pady=4)
        r += 1
        for variable in (self.v_ram, self.v_figh, self.v_dpi):
            variable.trace_add('write', lambda *_: self._update_allocation())
        self._update_allocation()
        self.after(2000, self._refresh_resources)
        return r

    @staticmethod
    def _ram_max(resources):
        return math.floor(max(0, resources.available-recruits.reserve_bytes(resources))
                          / recruits.GIB*4)/4

    def _refresh_resources(self):
        self._update_allocation()
        self.after(2000, self._refresh_resources)

    def _profile_key(self):
        paths = ((self.v_dir_a.get(), self.v_dir_b.get()) if self._comparison
                 else (self.v_geom.get(), self.v_data.get()))
        scalars = tuple(k for k, v in self.scalar_vars.items() if v.get())
        axes = tuple(k for k, v in self.axis_vars.items() if v.get())
        stamps = tuple(os.path.getmtime(p) if os.path.isfile(p) else None for p in paths)
        return (paths, scalars, axes, stamps, repr(cfg.CAR_BOUNDS_2D), repr(cfg.SCALARS),
                self.v_figh.get(), self.v_dpi.get())

    def _estimate_memory(self, bounds, scalar_count):
        if self._memory_profile is None or self._memory_profile[0] != self._profile_key():
            raise ValueError('Analyze RAM from input files to calculate RAM per recruit and recruit count.')
        profile = self._memory_profile[1]
        if 'estimate_bytes' in profile:
            return profile['estimate_bytes']
        return recruits.estimate_recruit_bytes(
            bounds, float(self.v_res.get()), scalar_count, **profile,
            fig_height=float(self.v_figh.get()), dpi=int(self.v_dpi.get()))

    def _update_allocation(self):
        if (not hasattr(self, 'allocation_lbl') or self._running or self._analyzing
                or self._updating_allocation):
            return
        self._updating_allocation = True
        try:
            resources = recruits.system_resources()
            ram_max = self._ram_max(resources)
            self.ram_slider.config(to=ram_max)
            if self.v_ram.get() > ram_max:
                self.v_ram.set(ram_max)
            if self.v_ram.get() <= 0:
                raise ValueError('No RAM budget selected or available after the system reserve.')
            bounds = [cfg.CAR_BOUNDS_2D[a] for a in range(3)
                      if self.axis_vars[cfg.AXIS_LABEL[a]].get()]
            scalar_count = sum(v.get() for v in self.scalar_vars.values())
            if not bounds or not scalar_count:
                raise ValueError('Select scalars and axes to estimate recruits.')
            estimate = self._estimate_memory(bounds, scalar_count)
            allocation = recruits.allocate_recruits(self.v_ram.get(), estimate,
                resources=resources, max_recruits=8 if self.v_parallel.get() else 1)
            text = (f'Estimated recruits: {allocation.recruits}  |  '
                    f'Budget: {allocation.budget_gib:.2f} GiB  |  '
                    f'~{allocation.per_recruit_gib:.2f} GiB per recruit\n'
                    f'Available: {resources.available/recruits.GIB:.1f} GiB  |  '
                    f'System reserve: {allocation.reserve_gib:.1f} GiB  |  '
                    f'CPU limit: {min(8, resources.cpus)} ({resources.cpus} logical CPUs available)\n'
                    f'Maximum budget now: {ram_max:.2f} GiB')
            if allocation.recruits == 0:
                text += '\nInsufficient budget or available RAM; increase budget or use a coarser resolution.'
            self.allocation_lbl.config(text=text, fg=ACCENT if allocation.recruits else '#C62828')
        except (ValueError, OverflowError, OSError, tk.TclError) as exc:
            self.allocation_lbl.config(text=str(exc), fg='#C62828')
        finally:
            self._updating_allocation = False

    def _analyze_ram(self):
        if self._running or self._analyzing:
            return
        key = self._profile_key()
        paths, scalars, axes = key[:3]
        if not all(paths) or not scalars or not axes:
            messagebox.showwarning('Missing inputs', 'Choose both input paths, scalars and axes.')
            return
        try:
            figh, dpi = float(self.v_figh.get()), int(self.v_dpi.get())
            if not math.isfinite(figh) or figh <= 0 or dpi <= 0:
                raise ValueError()
        except ValueError:
            messagebox.showerror('Invalid size', 'Use positive figure height and DPI.')
            return
        self._analyzing = True
        self.analyze_button.config(state='disabled')
        self.run_button.config(state='disabled')
        self.allocation_lbl.config(text='Analyzing input sizes...', fg=ACCENT)
        comparison = self._comparison
        def analyze():
            try:
                if comparison:
                    import compare
                    _, _, estimate = compare.prepare_comparison(*paths, scalars, axes, figh, dpi)
                    profile = dict(estimate_bytes=estimate)
                else:
                    import pipeline
                    profile = pipeline.inspect_recruit_memory(*paths, scalars, axes)
                self._analysis_events.put((key, profile, None))
            except Exception as exc:
                self._analysis_events.put((key, None, str(exc)))
        threading.Thread(target=analyze, daemon=True).start()
        self.after(100, self._poll_analysis)

    def _poll_analysis(self):
        try:
            key, profile, error = self._analysis_events.get_nowait()
        except queue.Empty:
            self.after(100, self._poll_analysis)
            return
        self._analyzing = False
        self.analyze_button.config(state='normal')
        self.run_button.config(state='normal')
        if error:
            self.allocation_lbl.config(text=f'RAM analysis failed: {error}', fg='#C62828')
        else:
            self._memory_profile = (key, profile)
            self._update_allocation()

    def _start_recruit_run(self):
        self._run_key = self._profile_key()
        self._running = True
        for widget in (self.run_button, self.analyze_button, self.ram_slider, self.parallel_check):
            widget.config(state='disabled')
        self.progress.config(value=0)
        self.status_lbl.config(text='Loading / staging...', fg=ACCENT)
        self.app.notebook.select(self.app.log_tab)
        self.after(100, self._poll_run)

    def _poll_run(self):
        while not self._events.empty():
            event = self._events.get_nowait()
            kind = event['event']
            if kind == 'allocation':
                profile = ({k: event[k] for k in ('point_count', 'geometry_bytes')}
                           if 'point_count' in event else dict(estimate_bytes=event['estimate_bytes']))
                self._memory_profile = (self._run_key, profile)
                self.allocation_lbl.config(text=f"Running with {event['recruits']} recruits  |  "
                    f"~{event['per_recruit_gib']:.2f} GiB per recruit  |  {event['total']} planes", fg=ACCENT)
                self.progress.config(maximum=event['total'], value=0)
            elif kind == 'plane':
                self.progress.config(value=event['processed'])
                self.status_lbl.config(text=f"{event['processed']}/{event['total']} planes | "
                                            f"{event['recruits']} recruits")
            elif kind in ('done', 'error'):
                self._running = False
                for widget in (self.run_button, self.analyze_button, self.ram_slider, self.parallel_check):
                    widget.config(state='normal')
                if kind == 'done':
                    self.status_lbl.config(text=f"Done: {event['completed']} saved, "
                                                f"{event['skipped']} skipped", fg='#2E7D32')
                else:
                    self.status_lbl.config(text='Error — see Log', fg='#C62828')
                    messagebox.showerror('Pipeline error', event['message'])
        if self._running:
            self.after(100, self._poll_run)
