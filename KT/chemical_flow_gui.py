"""
Chemical Process Flow Network - GUI Version

Interactive visual editor for creating chemical process flow networks.
- Draw squares (processes) and arrows (tubes)
- Define chemical components and track moles in each tube
- Set parameters (alpha per component, beta per tube/component)
- Alpha-eff applies to the overall process (from input to output)
- Solve and visualize moles of each component in all tubes
"""

import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from enum import Enum
import json


class ProcessType(Enum):
    INPUT = "input"
    OUTPUT = "output"
    REACTOR = "reactor"
    MIXER = "mixer"          # Multiple inputs → one output (combines flows)
    SEPARATOR = "separator"  # One input → multiple outputs (splits flow by beta)


# Colors for different process types
PROCESS_COLORS = {
    ProcessType.INPUT: "#90EE90",      # Light green
    ProcessType.OUTPUT: "#FFB6C1",     # Light pink
    ProcessType.REACTOR: "#87CEEB",    # Sky blue
    ProcessType.MIXER: "#F0E68C",      # Khaki
    ProcessType.SEPARATOR: "#DDA0DD",   # Plum
}


@dataclass
class Component:
    """Chemical component/species"""
    id: str
    name: str


@dataclass
class Process:
    id: str
    name: str
    process_type: ProcessType
    x: int
    y: int
    width: int = 100
    height: int = 60
    # Stoichiometric coefficients for reactors: component_id -> coefficient
    # Negative = reactant (consumed), Positive = product (produced)
    # e.g., {"A": -2, "B": -1, "C": 1, "D": 3} for 2A + B -> C + 3D
    stoich: Dict[str, float] = field(default_factory=dict)
    # Conversion (alpha): fraction of limiting reactant converted (0-1)
    # extent = conversion * (limiting reactant input / |stoich coeff|)
    conversion: float = 0.0
    # Beta per tube per component (for separators)
    beta_components: Dict[str, Dict[str, float]] = field(default_factory=dict)


@dataclass
class Tube:
    id: str
    from_process: str
    to_process: str
    beta: float = 1.0  # Split fraction (only used when alpha_eff = 0)
    moles: Dict[str, float] = field(default_factory=dict)  # Component -> moles
    
    def total_moles(self) -> float:
        return sum(self.moles.values())


class FlowNetworkGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Chemical Process Flow Network - Molar Flow Calculator")
        self.root.geometry("1400x900")
        
        # Data structures
        self.components: Dict[str, Component] = {}  # Chemical species
        self.processes: Dict[str, Process] = {}
        self.tubes: Dict[str, Tube] = {}
        self.input_moles: Dict[str, Dict[str, float]] = {}  # process_id -> comp_id -> moles
        self.alpha_eff: float = 0.0  # Overall process efficiency (0 = disabled)
        
        # GUI state
        self.selected_process = None
        self.drawing_tube = False
        self.tube_start_process = None
        self.process_counter = 0
        self.tube_counter = 0
        self.component_counter = 0
        
        # Undo history
        self.undo_stack = []
        self.max_undo = 50
        
        self.setup_gui()
        
    def setup_gui(self):
        """Set up the GUI layout"""
        # Main frame with horizontal paned window for resizable panels
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Horizontal PanedWindow for left panel | canvas | right panel
        self.main_paned = tk.PanedWindow(main_frame, orient=tk.HORIZONTAL, sashwidth=6, sashrelief=tk.RAISED)
        self.main_paned.pack(fill=tk.BOTH, expand=True)
        
        # Left panel - Tools (collapsible)
        self.left_panel = ttk.Frame(self.main_paned, width=180)
        self.left_panel_content = ttk.Frame(self.left_panel)
        self.left_panel_content.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Process buttons
        ttk.Label(self.left_panel_content, text="Add Process:", font=('Arial', 10, 'bold')).pack(pady=5)
        
        for ptype in ProcessType:
            btn = ttk.Button(self.left_panel_content, text=ptype.value.capitalize(),
                           command=lambda t=ptype: self.add_process_mode(t))
            btn.pack(fill=tk.X, pady=2)
        
        ttk.Separator(self.left_panel_content, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # Tube button
        ttk.Label(self.left_panel_content, text="Connections:", font=('Arial', 10, 'bold')).pack(pady=5)
        self.tube_btn = ttk.Button(self.left_panel_content, text="Draw Tube (Arrow)",
                                   command=self.toggle_tube_mode)
        self.tube_btn.pack(fill=tk.X, pady=2)
        
        ttk.Separator(self.left_panel_content, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # Actions
        ttk.Label(self.left_panel_content, text="Actions:", font=('Arial', 10, 'bold')).pack(pady=5)
        ttk.Button(self.left_panel_content, text="Manage Components", 
                  command=self.manage_components).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel_content, text="Set α_eff (Overall)", 
                  command=self.set_alpha_eff).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel_content, text="Solve Network", 
                  command=self.solve_network).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel_content, text="Clear All", 
                  command=self.clear_all).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel_content, text="Show Equations", 
                  command=self.show_equations).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel_content, text="Show Symbolic Eq.", 
                  command=self.show_symbolic_equations).pack(fill=tk.X, pady=2)
        
        ttk.Separator(self.left_panel_content, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # File operations
        ttk.Label(self.left_panel_content, text="File:", font=('Arial', 10, 'bold')).pack(pady=5)
        ttk.Button(self.left_panel_content, text="Save Network", 
                  command=self.save_network).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel_content, text="Load Network", 
                  command=self.load_network).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel_content, text="Load Example", 
                  command=self.load_example).pack(fill=tk.X, pady=2)
        
        self.main_paned.add(self.left_panel, minsize=20, width=180)
        
        # Canvas for drawing (center) with toggle buttons
        canvas_frame = ttk.Frame(self.main_paned)
        
        # Top bar with toggle buttons
        canvas_top_bar = ttk.Frame(canvas_frame)
        canvas_top_bar.pack(fill=tk.X)
        
        # Left panel toggle button
        self.left_visible = True
        self.toggle_left_btn = ttk.Button(canvas_top_bar, text="◀ Tools", width=10,
                                          command=self.toggle_left_panel)
        self.toggle_left_btn.pack(side=tk.LEFT, padx=2, pady=2)
        
        # Right panel toggle button
        self.right_visible = True
        self.toggle_right_btn = ttk.Button(canvas_top_bar, text="Results ▶", width=10,
                                           command=self.toggle_right_panel)
        self.toggle_right_btn.pack(side=tk.RIGHT, padx=2, pady=2)
        
        self.canvas = tk.Canvas(canvas_frame, bg='white', width=800, height=600)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        
        # Bind canvas events
        self.canvas.bind('<Button-1>', self.on_canvas_click)
        self.canvas.bind('<Button-3>', self.on_right_click)
        self.canvas.bind('<B1-Motion>', self.on_drag)
        self.canvas.bind('<ButtonRelease-1>', self.on_release)
        
        # Bind Ctrl+Z for undo
        self.root.bind('<Control-z>', self.undo)
        self.root.bind('<Control-Z>', self.undo)
        
        self.main_paned.add(canvas_frame, minsize=200)
        
        # Right panel - Results with collapsible sections
        self.right_panel = ttk.Frame(self.main_paned, width=280)
        
        # Results header
        results_header = ttk.Frame(self.right_panel)
        results_header.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(results_header, text="Results:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT)
        
        # Scrollable frame for results
        results_canvas = tk.Canvas(self.right_panel, highlightthickness=0)
        results_scrollbar = ttk.Scrollbar(self.right_panel, orient=tk.VERTICAL, command=results_canvas.yview)
        self.results_scroll_frame = ttk.Frame(results_canvas)
        
        self.results_scroll_frame.bind(
            "<Configure>",
            lambda e: results_canvas.configure(scrollregion=results_canvas.bbox("all"))
        )
        
        results_canvas.create_window((0, 0), window=self.results_scroll_frame, anchor="nw")
        results_canvas.configure(yscrollcommand=results_scrollbar.set)
        
        # Enable mouse wheel scrolling
        def on_mousewheel(event):
            results_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        results_canvas.bind_all("<MouseWheel>", on_mousewheel)
        
        results_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        results_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)
        
        self.main_paned.add(self.right_panel, minsize=20, width=280)
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready. Add processes and tubes to create network.")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.current_mode = None
        self.drag_process = None
        
        # Store reference to results canvas for scrolling
        self.results_canvas = results_canvas
        
        # Dictionary to store collapsible section states
        self.collapsed_sections = {}
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready. Add processes and tubes to create network.")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.current_mode = None
        self.drag_process = None
        
    def add_process_mode(self, process_type: ProcessType):
        """Enter mode to add a new process"""
        self.current_mode = ('add_process', process_type)
        self.drawing_tube = False
        self.status_var.set(f"Click on canvas to place {process_type.value}")
    
    def create_collapsible_section(self, parent, title, content_func, section_id, initially_open=False):
        """Create a collapsible section with header and content"""
        frame = ttk.Frame(parent)
        frame.pack(fill=tk.X, padx=2, pady=2)
        
        # Header frame with toggle button
        header = ttk.Frame(frame)
        header.pack(fill=tk.X)
        
        # Initialize collapse state
        if section_id not in self.collapsed_sections:
            self.collapsed_sections[section_id] = not initially_open
        
        is_collapsed = self.collapsed_sections[section_id]
        
        toggle_btn = ttk.Button(header, text="▶" if is_collapsed else "▼", width=2)
        toggle_btn.pack(side=tk.LEFT)
        
        title_label = ttk.Label(header, text=title, font=('Arial', 9, 'bold'))
        title_label.pack(side=tk.LEFT, padx=5)
        
        # Content frame
        content_frame = ttk.Frame(frame)
        if not is_collapsed:
            content_frame.pack(fill=tk.X, padx=15)
            content_func(content_frame)
        
        def toggle():
            self.collapsed_sections[section_id] = not self.collapsed_sections[section_id]
            if self.collapsed_sections[section_id]:
                content_frame.pack_forget()
                toggle_btn.config(text="▶")
            else:
                content_frame.pack(fill=tk.X, padx=15)
                # Clear and repopulate content
                for widget in content_frame.winfo_children():
                    widget.destroy()
                content_func(content_frame)
                toggle_btn.config(text="▼")
        
        toggle_btn.config(command=toggle)
        title_label.bind('<Button-1>', lambda e: toggle())
        
        return frame
    
    def toggle_left_panel(self):
        """Toggle visibility of left tools panel"""
        if self.left_visible:
            self.main_paned.forget(self.left_panel)
            self.toggle_left_btn.config(text="Tools ▶")
            self.left_visible = False
        else:
            # Re-add at position 0 (leftmost)
            self.main_paned.add(self.left_panel, before=self.main_paned.panes()[0], minsize=20, width=180)
            self.toggle_left_btn.config(text="◀ Tools")
            self.left_visible = True
    
    def toggle_right_panel(self):
        """Toggle visibility of right results panel"""
        if self.right_visible:
            self.main_paned.forget(self.right_panel)
            self.toggle_right_btn.config(text="◀ Results")
            self.right_visible = False
        else:
            # Re-add at the end (rightmost)
            self.main_paned.add(self.right_panel, minsize=20, width=280)
            self.toggle_right_btn.config(text="Results ▶")
            self.right_visible = True
        
    def toggle_tube_mode(self):
        """Toggle tube drawing mode"""
        self.drawing_tube = not self.drawing_tube
        self.current_mode = None
        if self.drawing_tube:
            self.status_var.set("Click on source process, then destination process")
            self.tube_btn.configure(text="Cancel Tube Drawing")
        else:
            self.status_var.set("Ready")
            self.tube_btn.configure(text="Draw Tube (Arrow)")
            self.tube_start_process = None
            
    def on_canvas_click(self, event):
        """Handle canvas click"""
        x, y = event.x, event.y
        
        # Check if clicked on a process
        clicked_process = self.get_process_at(x, y)
        
        if self.current_mode and self.current_mode[0] == 'add_process':
            # Adding new process
            process_type = self.current_mode[1]
            self.create_process(x, y, process_type)
            self.current_mode = None
            self.status_var.set("Ready")
            
        elif self.drawing_tube:
            if clicked_process:
                if self.tube_start_process is None:
                    self.tube_start_process = clicked_process
                    self.status_var.set(f"Now click on destination process")
                else:
                    if clicked_process != self.tube_start_process:
                        self.create_tube(self.tube_start_process, clicked_process)
                    self.tube_start_process = None
                    self.status_var.set("Click on source process, then destination process")
                    
        elif clicked_process:
            self.selected_process = clicked_process
            self.drag_process = clicked_process
            self.show_properties(clicked_process)
            
    def on_right_click(self, event):
        """Handle right click for context menu"""
        x, y = event.x, event.y
        clicked_process = self.get_process_at(x, y)
        clicked_tube = self.get_tube_at(x, y)
        
        menu = tk.Menu(self.root, tearoff=0)
        
        if clicked_process:
            menu.add_command(label=f"Edit {clicked_process}",
                           command=lambda: self.edit_process(clicked_process))
            menu.add_command(label=f"Delete {clicked_process}",
                           command=lambda: self.delete_process(clicked_process))
        elif clicked_tube:
            menu.add_command(label=f"Edit {clicked_tube}",
                           command=lambda: self.edit_tube(clicked_tube))
            menu.add_command(label=f"Delete {clicked_tube}",
                           command=lambda: self.delete_tube(clicked_tube))
            
        menu.post(event.x_root, event.y_root)
        
    def on_drag(self, event):
        """Handle dragging processes"""
        if self.drag_process and not self.drawing_tube:
            proc = self.processes.get(self.drag_process)
            if proc:
                proc.x = event.x
                proc.y = event.y
                self.redraw()
                
    def on_release(self, event):
        """Handle mouse release"""
        self.drag_process = None
        
    def get_process_at(self, x: int, y: int) -> Optional[str]:
        """Get process ID at coordinates"""
        for pid, proc in self.processes.items():
            if (proc.x - proc.width//2 <= x <= proc.x + proc.width//2 and
                proc.y - proc.height//2 <= y <= proc.y + proc.height//2):
                return pid
        return None
    
    def get_tube_at(self, x: int, y: int) -> Optional[str]:
        """Get tube ID near coordinates (works with routed paths)"""
        for tid, tube in self.tubes.items():
            from_proc = self.processes.get(tube.from_process)
            to_proc = self.processes.get(tube.to_process)
            if from_proc and to_proc:
                # Get the route points for this tube
                points = self.calculate_tube_route(tube, from_proc, to_proc)
                
                # Check distance to each segment
                for i in range(len(points) - 1):
                    x1, y1 = points[i]
                    x2, y2 = points[i + 1]
                    dist = self.point_to_line_distance(x, y, x1, y1, x2, y2)
                    if dist < 10:
                        return tid
        return None
        
    def point_to_line_distance(self, px, py, x1, y1, x2, y2) -> float:
        """Calculate distance from point to line segment"""
        dx = x2 - x1
        dy = y2 - y1
        
        if dx == 0 and dy == 0:
            return np.sqrt((px - x1)**2 + (py - y1)**2)
            
        t = max(0, min(1, ((px - x1)*dx + (py - y1)*dy) / (dx*dx + dy*dy)))
        
        proj_x = x1 + t * dx
        proj_y = y1 + t * dy
        
        return np.sqrt((px - proj_x)**2 + (py - proj_y)**2)
        
    def create_process(self, x: int, y: int, process_type: ProcessType):
        """Create a new process"""
        self.save_undo_state()
        self.process_counter += 1
        pid = f"P{self.process_counter}"
        
        name = simpledialog.askstring("Process Name", f"Enter name for {process_type.value}:",
                                     initialvalue=f"{process_type.value.capitalize()} {self.process_counter}",
                                     parent=self.root)
        if not name:
            name = f"{process_type.value.capitalize()} {self.process_counter}"
        
        stoich = {}
        conversion = 0.0
        if process_type == ProcessType.REACTOR:
            # Ask for stoichiometric coefficients
            if self.components:
                messagebox.showinfo("Stoichiometry", 
                    "Enter stoichiometric coefficients:\n" +
                    "• Negative for reactants (consumed)\n" +
                    "• Positive for products (produced)\n" +
                    "• Zero if not involved\n\n" +
                    "Example: 2A + B → C + 3D\n" +
                    "A=-2, B=-1, C=+1, D=+3",
                    parent=self.root)
                for comp_id, comp in self.components.items():
                    coeff = simpledialog.askfloat(f"Coefficient for {comp.name}", 
                                             f"Stoichiometric coeff for {comp.name}:\n(negative=reactant, positive=product, 0=not involved)",
                                             initialvalue=0.0, parent=self.root)
                    if coeff is not None:
                        stoich[comp_id] = coeff
                
                # Ask for conversion (alpha)
                conversion = simpledialog.askfloat("Conversion (α)", 
                    "Enter conversion (0-1):\n" +
                    "Fraction of limiting reactant that reacts.\n" +
                    "(e.g., 0.8 = 80% conversion)",
                    initialvalue=0.8, minvalue=0.0, maxvalue=1.0, parent=self.root)
                if conversion is None:
                    conversion = 0.8
                
        if process_type == ProcessType.INPUT:
            # Ask for input moles per component
            if self.components:
                self.input_moles[pid] = {}
                for comp_id, comp in self.components.items():
                    moles = simpledialog.askfloat(f"Input {comp.name}", 
                                                 f"Enter input moles of {comp.name}:",
                                                 initialvalue=0.0, minvalue=0.0,
                                                 parent=self.root)
                    if moles and moles > 0:
                        self.input_moles[pid][comp_id] = moles
            else:
                messagebox.showinfo("Info", "Add components first using 'Manage Components'", parent=self.root)
        
        self.processes[pid] = Process(pid, name, process_type, x, y, stoich=stoich, conversion=conversion)
        self.redraw()
        self.status_var.set(f"Created {name}")
        
    def create_tube(self, from_id: str, to_id: str):
        """Create a new tube"""
        self.save_undo_state()
        self.tube_counter += 1
        tid = f"T{self.tube_counter}"
        
        # Check if this is a recirculation (going back to an earlier process)
        # Default beta = 1.0 (can be changed by editing the tube)
        self.tubes[tid] = Tube(tid, from_id, to_id, beta=1.0)
        
        self.redraw()
        self.status_var.set(f"Created tube from {from_id} to {to_id} (β=1.0)")
        
    def edit_process(self, pid: str):
        """Edit process properties"""
        proc = self.processes.get(pid)
        if not proc:
            return
            
        # Save state for undo
        self.save_undo_state()
        
        # Create edit dialog
        dialog = tk.Toplevel(self.root)
        dialog.title(f"Edit {pid} - {proc.name}")
        dialog.geometry("400x450")
        dialog.transient(self.root)  # Keep on top of main window
        dialog.grab_set()  # Modal dialog
        
        row = 0
        ttk.Label(dialog, text="Name:").grid(row=row, column=0, padx=5, pady=5, sticky='e')
        name_var = tk.StringVar(value=proc.name)
        first_entry = ttk.Entry(dialog, textvariable=name_var)
        first_entry.grid(row=row, column=1, padx=5, pady=5)
        # Auto-focus on first entry
        dialog.after(100, lambda: first_entry.focus_set())
        dialog.after(100, lambda: first_entry.selection_range(0, tk.END))
        row += 1
        
        # Stoichiometric coefficients for reactors
        stoich_vars = {}
        conversion_var = None
        if proc.process_type == ProcessType.REACTOR and self.components:
            ttk.Label(dialog, text="Stoichiometric coefficients:", font=('Arial', 9, 'bold')).grid(
                row=row, column=0, columnspan=2, pady=5)
            row += 1
            ttk.Label(dialog, text="(negative=reactant, positive=product)", font=('Arial', 8)).grid(
                row=row, column=0, columnspan=2)
            row += 1
            for comp_id, comp in self.components.items():
                ttk.Label(dialog, text=f"  {comp.name}:").grid(row=row, column=0, padx=5, pady=2, sticky='e')
                stoich_vars[comp_id] = tk.DoubleVar(value=proc.stoich.get(comp_id, 0.0))
                ttk.Entry(dialog, textvariable=stoich_vars[comp_id], width=10).grid(row=row, column=1, padx=5, pady=2, sticky='w')
                row += 1
            
            # Conversion (alpha)
            ttk.Separator(dialog, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky='ew', pady=5)
            row += 1
            ttk.Label(dialog, text="Conversion (α):", font=('Arial', 9, 'bold')).grid(
                row=row, column=0, padx=5, pady=2, sticky='e')
            conversion_var = tk.DoubleVar(value=proc.conversion)
            ttk.Entry(dialog, textvariable=conversion_var, width=10).grid(row=row, column=1, padx=5, pady=2, sticky='w')
            row += 1
            ttk.Label(dialog, text="(0-1, fraction of limiting reactant)", font=('Arial', 8)).grid(
                row=row, column=0, columnspan=2)
            row += 1
        
        # Input moles for input processes
        input_vars = {}
        if proc.process_type == ProcessType.INPUT and self.components:
            ttk.Label(dialog, text="Input moles per component:", font=('Arial', 9, 'bold')).grid(
                row=row, column=0, columnspan=2, pady=5)
            row += 1
            for comp_id, comp in self.components.items():
                ttk.Label(dialog, text=f"  {comp.name}:").grid(row=row, column=0, padx=5, pady=2, sticky='e')
                current = self.input_moles.get(pid, {}).get(comp_id, 0.0)
                input_vars[comp_id] = tk.DoubleVar(value=current)
                ttk.Entry(dialog, textvariable=input_vars[comp_id], width=10).grid(row=row, column=1, padx=5, pady=2, sticky='w')
                row += 1
        
        def save():
            proc.name = name_var.get()
            # Save stoichiometry
            for comp_id, var in stoich_vars.items():
                proc.stoich[comp_id] = var.get()
            # Save conversion for reactors
            if conversion_var is not None:
                proc.conversion = conversion_var.get()
            # Save input moles
            if proc.process_type == ProcessType.INPUT:
                if pid not in self.input_moles:
                    self.input_moles[pid] = {}
                for comp_id, var in input_vars.items():
                    if var.get() > 0:
                        self.input_moles[pid][comp_id] = var.get()
            self.redraw()
            dialog.destroy()
            
        ttk.Button(dialog, text="Save", command=save).grid(row=row, column=0, columnspan=2, pady=10)
        
    def edit_tube(self, tid: str):
        """Edit tube properties"""
        tube = self.tubes.get(tid)
        if not tube:
            return
        
        # Save state for undo
        self.save_undo_state()
        
        # Create edit dialog
        dialog = tk.Toplevel(self.root)
        dialog.title(f"Edit Tube {tid}")
        dialog.geometry("400x400")
        dialog.transient(self.root)
        dialog.grab_set()
        
        row = 0
        
        ttk.Label(dialog, text="Beta (split fraction 0-1):").grid(row=row, column=0, padx=5, pady=5, sticky='e')
        beta_var = tk.DoubleVar(value=tube.beta)
        first_entry = ttk.Entry(dialog, textvariable=beta_var, width=10)
        first_entry.grid(row=row, column=1, padx=5, pady=5, sticky='w')
        # Auto-focus on first entry
        dialog.after(100, lambda: first_entry.focus_set())
        dialog.after(100, lambda: first_entry.selection_range(0, tk.END))
        row += 1
        
        ttk.Label(dialog, text="(Only used when α_eff = 0)", font=('Arial', 8)).grid(
            row=row, column=0, columnspan=2)
        row += 1
        
        # Component-specific betas for separator outputs
        from_process = self.processes.get(tube.from_process)
        beta_comp_vars = {}
        if from_process and from_process.process_type == ProcessType.SEPARATOR and self.components:
            ttk.Label(dialog, text="Beta per component (0 = unknown):", font=('Arial', 9, 'bold')).grid(
                row=row, column=0, columnspan=2, pady=5)
            row += 1
            ttk.Label(dialog, text="(only used when α_eff = 0)", font=('Arial', 8)).grid(
                row=row, column=0, columnspan=2)
            row += 1
            for comp_id, comp in self.components.items():
                ttk.Label(dialog, text=f"  β({comp.name}):").grid(row=row, column=0, padx=5, pady=2, sticky='e')
                # Default to 0 (unknown) if not previously set, not tube.beta
                current = from_process.beta_components.get(tid, {}).get(comp_id, 0.0)
                beta_comp_vars[comp_id] = tk.DoubleVar(value=current)
                ttk.Entry(dialog, textvariable=beta_comp_vars[comp_id], width=10).grid(row=row, column=1, padx=5, pady=2, sticky='w')
                row += 1
        
        def save():
            tube.beta = beta_var.get()
            if from_process and from_process.process_type == ProcessType.SEPARATOR:
                if tid not in from_process.beta_components:
                    from_process.beta_components[tid] = {}
                for comp_id, var in beta_comp_vars.items():
                    from_process.beta_components[tid][comp_id] = var.get()
            self.redraw()
            dialog.destroy()
            
        ttk.Button(dialog, text="Save", command=save).grid(row=row, column=0, columnspan=2, pady=10)
    
    def save_undo_state(self):
        """Save current state to undo stack"""
        import copy
        state = {
            'components': copy.deepcopy({k: {'id': v.id, 'name': v.name} for k, v in self.components.items()}),
            'processes': copy.deepcopy({k: {
                'id': v.id, 'name': v.name, 'type': v.process_type.value,
                'x': v.x, 'y': v.y, 'stoich': v.stoich, 
                'conversion': v.conversion, 'beta_components': v.beta_components
            } for k, v in self.processes.items()}),
            'tubes': copy.deepcopy({k: {
                'id': v.id, 'from': v.from_process, 'to': v.to_process, 'beta': v.beta
            } for k, v in self.tubes.items()}),
            'input_moles': copy.deepcopy(self.input_moles),
            'alpha_eff': self.alpha_eff,
            'process_counter': self.process_counter,
            'tube_counter': self.tube_counter,
            'component_counter': self.component_counter
        }
        self.undo_stack.append(state)
        # Limit stack size
        if len(self.undo_stack) > self.max_undo:
            self.undo_stack.pop(0)
    
    def undo(self, event=None):
        """Undo last action"""
        if not self.undo_stack:
            self.status_var.set("Nothing to undo")
            return
        
        state = self.undo_stack.pop()
        
        # Restore components
        self.components.clear()
        for cid, cdata in state['components'].items():
            self.components[cid] = Component(cdata['id'], cdata['name'])
        
        # Restore processes
        self.processes.clear()
        for pid, pdata in state['processes'].items():
            self.processes[pid] = Process(
                pdata['id'], pdata['name'], ProcessType(pdata['type']),
                pdata['x'], pdata['y'],
                stoich=pdata.get('stoich', {}),
                conversion=pdata.get('conversion', 0.0),
                beta_components=pdata.get('beta_components', {})
            )
        
        # Restore tubes
        self.tubes.clear()
        for tid, tdata in state['tubes'].items():
            self.tubes[tid] = Tube(
                tdata['id'], tdata['from'], tdata['to'], tdata.get('beta', 1.0)
            )
        
        self.input_moles = state['input_moles']
        self.alpha_eff = state['alpha_eff']
        self.process_counter = state['process_counter']
        self.tube_counter = state['tube_counter']
        self.component_counter = state['component_counter']
        
        self.redraw()
        self.status_var.set("Undo successful")
            
    def delete_process(self, pid: str):
        """Delete a process and its connected tubes"""
        self.save_undo_state()
        if pid in self.processes:
            del self.processes[pid]
            if pid in self.input_moles:
                del self.input_moles[pid]
            # Delete connected tubes
            to_delete = [tid for tid, tube in self.tubes.items() 
                        if tube.from_process == pid or tube.to_process == pid]
            for tid in to_delete:
                del self.tubes[tid]
            self.redraw()
            
    def delete_tube(self, tid: str):
        """Delete a tube"""
        self.save_undo_state()
        if tid in self.tubes:
            del self.tubes[tid]
            self.redraw()
            
    def show_properties(self, pid: str):
        """Show properties of selected process"""
        proc = self.processes.get(pid)
        if not proc:
            return
            
        # Clear properties frame
        for widget in self.props_frame.winfo_children():
            widget.destroy()
            
        ttk.Label(self.props_frame, text=f"ID: {pid}").pack(anchor=tk.W)
        ttk.Label(self.props_frame, text=f"Name: {proc.name}").pack(anchor=tk.W)
        ttk.Label(self.props_frame, text=f"Type: {proc.process_type.value}").pack(anchor=tk.W)
        
        if proc.process_type == ProcessType.REACTOR and proc.stoich:
            ttk.Label(self.props_frame, text="Stoichiometry:").pack(anchor=tk.W)
            for comp_id, coeff in proc.stoich.items():
                sign = "+" if coeff > 0 else ""
                ttk.Label(self.props_frame, text=f"  {comp_id}: {sign}{coeff:.2f}").pack(anchor=tk.W)
        
        if pid in self.input_moles:
            ttk.Label(self.props_frame, text="Input moles:").pack(anchor=tk.W)
            for comp_id, moles in self.input_moles[pid].items():
                ttk.Label(self.props_frame, text=f"  {comp_id}: {moles:.2f}").pack(anchor=tk.W)
    
    def manage_components(self):
        """Dialog to add/edit chemical components"""
        self.save_undo_state()
        
        dialog = tk.Toplevel(self.root)
        dialog.title("Manage Chemical Components")
        dialog.geometry("400x350")
        dialog.transient(self.root)
        
        # Listbox for components
        ttk.Label(dialog, text="Components:", font=('Arial', 10, 'bold')).pack(pady=5)
        
        list_frame = ttk.Frame(dialog)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        listbox = tk.Listbox(list_frame, height=10)
        listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=listbox.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        listbox.configure(yscrollcommand=scrollbar.set)
        
        def refresh_list():
            listbox.delete(0, tk.END)
            for cid, comp in self.components.items():
                listbox.insert(tk.END, cid)  # Just show the name/ID
        
        refresh_list()
        
        # Add component frame
        add_frame = ttk.LabelFrame(dialog, text="Add New Component")
        add_frame.pack(fill=tk.X, padx=10, pady=5)
        
        inner_frame = ttk.Frame(add_frame)
        inner_frame.pack(padx=5, pady=5)
        
        ttk.Label(inner_frame, text="Name:").grid(row=0, column=0, padx=2, sticky='e')
        name_var = tk.StringVar()
        name_entry = ttk.Entry(inner_frame, textvariable=name_var, width=20)
        name_entry.grid(row=0, column=1, padx=2)
        # Auto-focus on entry
        dialog.after(100, lambda: name_entry.focus_set())
        
        def add_component():
            name = name_var.get().strip()
            if name:
                # Use the name directly as the ID (e.g., H2, N2, NH3)
                cid = name
                if cid in self.components:
                    messagebox.showwarning("Warning", f"Component '{cid}' already exists!", parent=dialog)
                    return
                    
                self.components[cid] = Component(cid, name)
                refresh_list()
                name_var.set("")
                self.status_var.set(f"Added component {cid}")
        
        def delete_selected():
            sel = listbox.curselection()
            if sel:
                cid = list(self.components.keys())[sel[0]]
                del self.components[cid]
                refresh_list()
                self.status_var.set(f"Deleted component {cid}")
        
        ttk.Button(inner_frame, text="Add", command=add_component).grid(row=0, column=2, padx=10)
        
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(pady=10)
        ttk.Button(btn_frame, text="Delete Selected", command=delete_selected).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Close", command=dialog.destroy).pack(side=tk.LEFT, padx=5)
    
    def set_alpha_eff(self):
        """Set overall process efficiency"""
        val = simpledialog.askfloat("Overall Process Efficiency", 
                                   "Enter α_eff (0-1):\\n" +
                                   "Fraction of limiting reactant converted overall.\\n" +
                                   "α_eff=1.0 → 100% converted (0 mol output)\\n" +
                                   "α_eff=0.5 → 50% converted\\n" +
                                   "α_eff=0 → disabled (no constraint)",
                                   initialvalue=self.alpha_eff, minvalue=0.0, maxvalue=1.0,
                                   parent=self.root)
        if val is not None:
            self.save_undo_state()
            self.alpha_eff = val
            if val > 0:
                self.status_var.set(f"Set overall α_eff = {val:.4f} (enabled)")
            else:
                self.status_var.set("α_eff disabled (set to 0)")
            
    def redraw(self):
        """Redraw the entire canvas"""
        self.canvas.delete("all")
        
        # Draw tubes first (behind processes)
        for tid, tube in self.tubes.items():
            self.draw_tube(tube)
            
        # Draw processes
        for pid, proc in self.processes.items():
            self.draw_process(proc)
            
    def draw_process(self, proc: Process):
        """Draw a process on the canvas with auto-sizing based on content"""
        x, y = proc.x, proc.y
        
        color = PROCESS_COLORS.get(proc.process_type, "#FFFFFF")
        
        # Calculate required width based on content
        name_width = len(proc.name) * 8 + 20  # Approximate width for name
        
        # For reactors, calculate reaction string width
        rxn_str = ""
        if proc.process_type == ProcessType.REACTOR and proc.stoich:
            reactants = [f"{abs(v):.0f}{k}" for k, v in proc.stoich.items() if v < 0]
            products = [f"{v:.0f}{k}" for k, v in proc.stoich.items() if v > 0]
            rxn_str = " + ".join(reactants) + " → " + " + ".join(products) if products else ""
            rxn_width = len(rxn_str) * 6 + 20  # Approximate width for reaction
        else:
            rxn_width = 0
        
        # Set width to fit content (minimum 80, max based on content)
        min_width = 80
        content_width = max(name_width, rxn_width)
        w = max(min_width, content_width)
        
        # Height: taller if showing reaction
        h = 70 if rxn_str else 50
        
        # Update process dimensions for hit detection
        proc.width = w
        proc.height = h
        
        # Draw rectangle
        self.canvas.create_rectangle(x - w//2, y - h//2, x + w//2, y + h//2,
                                    fill=color, outline="black", width=2)
        
        # Draw name
        self.canvas.create_text(x, y - 10 if rxn_str else y, text=proc.name, font=('Arial', 9, 'bold'))
        
        # Draw stoichiometry for reactors (full string, no truncation)
        if rxn_str:
            self.canvas.create_text(x, y + 12, text=rxn_str, font=('Arial', 8))
        
        # Draw ID
        self.canvas.create_text(x - w//2 + 15, y - h//2 + 10, text=proc.id, 
                               font=('Arial', 8), fill="gray")
                               
    def draw_tube(self, tube: Tube):
        """Draw a tube (arrow) on the canvas with smart routing"""
        from_proc = self.processes.get(tube.from_process)
        to_proc = self.processes.get(tube.to_process)
        
        if not from_proc or not to_proc:
            return
        
        # Get the route points for this tube
        points = self.calculate_tube_route(tube, from_proc, to_proc)
        
        tube_color = "blue"
        
        # Draw the polyline with arrow at the end
        if len(points) >= 2:
            # Draw line segments
            for i in range(len(points) - 1):
                x1, y1 = points[i]
                x2, y2 = points[i + 1]
                is_last = (i == len(points) - 2)
                self.canvas.create_line(x1, y1, x2, y2, 
                                        arrow=tk.LAST if is_last else None, 
                                        width=2, fill=tube_color)
            
            # Draw label at midpoint of route
            mid_idx = len(points) // 2
            if len(points) % 2 == 0:
                mx = (points[mid_idx-1][0] + points[mid_idx][0]) / 2
                my = (points[mid_idx-1][1] + points[mid_idx][1]) / 2
            else:
                mx, my = points[mid_idx]
            
            label = f"{tube.id}\nβ={tube.beta:.2f}"
            if tube.moles:
                total = tube.total_moles()
                label += f"\n{total:.2f} mol"
            self.canvas.create_text(mx, my - 15, text=label, font=('Arial', 8), fill=tube_color)
    
    def calculate_tube_route(self, tube: Tube, from_proc: Process, to_proc: Process) -> List[Tuple[int, int]]:
        """Calculate the route points for a tube, avoiding other processes"""
        
        # Check if from a separator going to a mixer - use top/bottom edge
        from_is_separator = from_proc.process_type == ProcessType.SEPARATOR
        to_is_mixer = to_proc.process_type == ProcessType.MIXER
        use_vertical_routing = from_is_separator and to_is_mixer
        
        # Check if this is a recirculation (going backward in x direction)
        is_recirculation = to_proc.x < from_proc.x
        
        # For recirculation, separator->mixer, or blocked paths, use rectangular routing
        if is_recirculation or use_vertical_routing or self.direct_path_blocked(from_proc, to_proc):
            return self.calculate_rectangular_route(tube, from_proc, to_proc, is_recirculation, use_vertical_routing, to_is_mixer)
        else:
            # Simple direct connection
            return self.calculate_direct_route(from_proc, to_proc)
    
    def calculate_direct_route(self, from_proc: Process, to_proc: Process) -> List[Tuple[int, int]]:
        """Calculate direct route between two processes"""
        x1, y1 = self.get_edge_point(from_proc, to_proc)
        x2, y2 = self.get_edge_point(to_proc, from_proc)
        return [(x1, y1), (x2, y2)]
    
    def calculate_rectangular_route(self, tube: Tube, from_proc: Process, to_proc: Process, 
                                    is_recirculation: bool, use_vertical_exit: bool = False,
                                    to_is_mixer: bool = False) -> List[Tuple[int, int]]:
        """Calculate a rectangular route that goes around other processes"""
        points = []
        
        # Determine routing direction based on relative positions
        margin = 40  # Distance to route around processes
        
        # Get all process bounding boxes for collision detection
        all_procs = list(self.processes.values())
        
        # Find the vertical range of processes between from and to
        min_y = float('inf')
        max_y = float('-inf')
        for p in all_procs:
            if min(from_proc.x, to_proc.x) - 50 <= p.x <= max(from_proc.x, to_proc.x) + 50:
                min_y = min(min_y, p.y - p.height//2)
                max_y = max(max_y, p.y + p.height//2)
        
        if min_y == float('inf'):
            min_y = min(from_proc.y, to_proc.y) - 50
            max_y = max(from_proc.y, to_proc.y) + 50
        
        # Decide to go above or below
        go_below = (from_proc.y + to_proc.y) / 2 < 300
        
        if is_recirculation:
            # Recirculation: exit from bottom/top, go around, enter from bottom/top
            if go_below:
                y_route = max_y + margin
                y_start = from_proc.y + from_proc.height//2  # Exit from bottom
                y_end = to_proc.y + to_proc.height//2        # Enter from bottom
            else:
                y_route = min_y - margin
                y_start = from_proc.y - from_proc.height//2  # Exit from top
                y_end = to_proc.y - to_proc.height//2        # Enter from top
            
            x_start = from_proc.x
            x_end = to_proc.x
            
            points = [
                (x_start, y_start),           # Exit from_proc (top/bottom)
                (x_start, y_route),           # Go up/down
                (x_end, y_route),             # Go horizontally
                (x_end, y_end)                # Enter to_proc (top/bottom)
            ]
        elif use_vertical_exit and to_is_mixer:
            # Separator -> Mixer: exit from top or bottom based on tube index
            outgoing_tubes = [t for t in self.tubes.values() if t.from_process == from_proc.id]
            tube_index = next((i for i, t in enumerate(outgoing_tubes) if t.id == tube.id), 0)
            
            # First output goes down, second goes up, etc.
            if tube_index % 2 == 0:
                y_start = from_proc.y + from_proc.height//2  # Bottom
                y_route = max_y + margin
                y_end = to_proc.y + to_proc.height//2        # Enter mixer from bottom
            else:
                y_start = from_proc.y - from_proc.height//2  # Top
                y_route = min_y - margin
                y_end = to_proc.y - to_proc.height//2        # Enter mixer from top
            
            x_start = from_proc.x
            x_end = to_proc.x
            
            points = [
                (x_start, y_start),           # Exit from separator (top/bottom)
                (x_start, y_route),           # Go up/down
                (x_end, y_route),             # Go horizontally
                (x_end, y_end)                # Enter mixer (top/bottom)
            ]
        else:
            # Forward direction but path is blocked - route around
            x1, y1 = self.get_edge_point(from_proc, to_proc)
            x2, y2 = self.get_edge_point(to_proc, from_proc)
            
            # Check if there's a process in the way
            blocking_proc = self.find_blocking_process(from_proc, to_proc)
            
            if blocking_proc:
                # Route around the blocking process
                bp = blocking_proc
                
                # Go above or below the blocking process
                if from_proc.y < bp.y:
                    y_route = bp.y - bp.height//2 - margin
                else:
                    # Go below
                    y_route = bp.y + bp.height//2 + margin
                
                points = [
                    (x1, y1),
                    (x1, y_route),
                    (x2, y_route),
                    (x2, y2)
                ]
            else:
                points = [(x1, y1), (x2, y2)]
        
        return points
    
    def direct_path_blocked(self, from_proc: Process, to_proc: Process) -> bool:
        """Check if a direct path between two processes is blocked by another process"""
        return self.find_blocking_process(from_proc, to_proc) is not None
    
    def find_blocking_process(self, from_proc: Process, to_proc: Process) -> Optional[Process]:
        """Find a process that blocks the direct path between from and to"""
        x1, y1 = from_proc.x, from_proc.y
        x2, y2 = to_proc.x, to_proc.y
        
        for proc in self.processes.values():
            if proc.id == from_proc.id or proc.id == to_proc.id:
                continue
            
            # Check if the line from (x1,y1) to (x2,y2) intersects the process rectangle
            if self.line_intersects_rect(x1, y1, x2, y2, 
                                         proc.x - proc.width//2 - 10,
                                         proc.y - proc.height//2 - 10,
                                         proc.x + proc.width//2 + 10,
                                         proc.y + proc.height//2 + 10):
                return proc
        return None
    
    def line_intersects_rect(self, x1, y1, x2, y2, rx1, ry1, rx2, ry2) -> bool:
        """Check if a line segment intersects a rectangle"""
        # Check if line segment intersects any of the 4 edges of the rectangle
        edges = [
            (rx1, ry1, rx2, ry1),  # Top
            (rx2, ry1, rx2, ry2),  # Right
            (rx1, ry2, rx2, ry2),  # Bottom
            (rx1, ry1, rx1, ry2),  # Left
        ]
        
        for ex1, ey1, ex2, ey2 in edges:
            if self.lines_intersect(x1, y1, x2, y2, ex1, ey1, ex2, ey2):
                return True
        
        # Also check if the line is completely inside the rectangle
        if rx1 <= x1 <= rx2 and ry1 <= y1 <= ry2:
            return True
        if rx1 <= x2 <= rx2 and ry1 <= y2 <= ry2:
            return True
            
        return False
    
    def lines_intersect(self, x1, y1, x2, y2, x3, y3, x4, y4) -> bool:
        """Check if two line segments intersect"""
        def ccw(ax, ay, bx, by, cx, cy):
            return (cy - ay) * (bx - ax) > (by - ay) * (cx - ax)
        
        return (ccw(x1, y1, x3, y3, x4, y4) != ccw(x2, y2, x3, y3, x4, y4) and
                ccw(x1, y1, x2, y2, x3, y3) != ccw(x1, y1, x2, y2, x4, y4))
        
    def get_edge_point(self, from_proc: Process, to_proc: Process) -> Tuple[int, int]:
        """Get the point on the edge of from_proc facing to_proc"""
        dx = to_proc.x - from_proc.x
        dy = to_proc.y - from_proc.y
        
        if dx == 0 and dy == 0:
            return from_proc.x, from_proc.y
            
        # Determine which edge to use
        w, h = from_proc.width // 2, from_proc.height // 2
        
        if abs(dx) * h > abs(dy) * w:
            # Horizontal edge
            if dx > 0:
                return from_proc.x + w, from_proc.y + int(dy * w / abs(dx))
            else:
                return from_proc.x - w, from_proc.y - int(dy * w / abs(dx))
        else:
            # Vertical edge
            if dy > 0:
                return from_proc.x + int(dx * h / abs(dy)), from_proc.y + h
            else:
                return from_proc.x - int(dx * h / abs(dy)), from_proc.y - h
                
    def solve_network(self):
        """Solve the flow network - calculates moles of each component in each tube"""
        if not self.processes or not self.tubes:
            messagebox.showwarning("Warning", "Please create a network first!", parent=self.root)
            return
        
        if not self.components:
            messagebox.showwarning("Warning", "Please define components first using 'Manage Components'!", parent=self.root)
            return
        
        # Validate separator betas sum to 1 (only when alpha_eff is NOT specified)
        if self.alpha_eff == 0:
            warnings = []
            for proc in self.processes.values():
                if proc.process_type == ProcessType.SEPARATOR:
                    outgoing = [t for t in self.tubes.values() if t.from_process == proc.id]
                    if len(outgoing) > 1:
                        for comp_id in self.components.keys():
                            beta_sum = 0.0
                            for tube in outgoing:
                                if tube.id in proc.beta_components:
                                    beta_sum += proc.beta_components[tube.id].get(comp_id, tube.beta)
                                else:
                                    beta_sum += tube.beta
                            if abs(beta_sum - 1.0) > 0.01:
                                warnings.append(f"  {proc.name}: β sum for {comp_id} = {beta_sum:.2f} (should be 1.0)")
        
            if warnings:
                msg = "WARNING: Separator betas don't sum to 1.0!\\nThis will cause mass imbalance:\\n" + "\\n".join(warnings)
                if not messagebox.askyesno("Beta Warning", msg + "\\n\\nContinue anyway?", parent=self.root):
                    return
            
        try:
            A, b, variables = self.build_equation_system()
            
            if len(A) == 0:
                messagebox.showwarning("Warning", "No equations to solve!", parent=self.root)
                return
                
            n_equations, n_variables = A.shape
            
            # Remove redundant equations
            A_reduced, b_reduced = self._remove_redundant_equations(A, b)
            
            # Solve
            if len(A_reduced) == n_variables:
                try:
                    solution = np.linalg.solve(A_reduced, b_reduced)
                except np.linalg.LinAlgError:
                    solution, _, _, _ = np.linalg.lstsq(A_reduced, b_reduced, rcond=None)
            else:
                solution, _, _, _ = np.linalg.lstsq(A_reduced, b_reduced, rcond=None)
            
            # Check for negative mole values and warn user
            negative_vars = []
            for i, var in enumerate(variables):
                if var[0] == 'tube' and solution[i] < -1e-6:
                    negative_vars.append((var[1], var[2], solution[i]))
            
            if negative_vars:
                # This indicates the system may be under-constrained or physically impossible
                msg = "Warning: Some intermediate flows have negative values:\\n"
                for tid, cid, val in negative_vars[:5]:  # Show first 5
                    msg += f"  {tid} {cid}: {val:.4f}\\n"
                if len(negative_vars) > 5:
                    msg += f"  ... and {len(negative_vars) - 5} more\\n"
                msg += "\\nThis may indicate:\\n"
                msg += "- The system is under-constrained (more unknowns than equations)\\n"
                msg += "- The α_eff target may not be physically achievable\\n"
                msg += "- The per-pass conversion α may be too low for the target α_eff"
                messagebox.showwarning("Negative Flows", msg, parent=self.root)
            
            # Clear old results
            for tube in self.tubes.values():
                tube.moles = {}
                
            # Store results (no multiplication - alpha_eff is handled in equations)
            for i, var in enumerate(variables):
                var_type = var[0]
                if var_type == 'tube':
                    _, tube_id, comp_id = var
                    self.tubes[tube_id].moles[comp_id] = solution[i]
                
            # Display results
            self.display_results()
            self.redraw()
            self.status_var.set("Network solved successfully!")
            
        except Exception as e:
            import traceback
            messagebox.showerror("Error", f"Failed to solve: {str(e)}\\n{traceback.format_exc()}", parent=self.root)
    
    def _remove_redundant_equations(self, A: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Remove linearly dependent equations"""
        if len(A) == 0:
            return A, b
        Q, R = np.linalg.qr(A.T)
        independent = np.abs(np.diag(R)) > 1e-10
        independent_indices = []
        rank = 0
        for i in range(len(A)):
            if rank < len(independent) and independent[rank]:
                independent_indices.append(i)
            rank += 1
            if rank >= len(independent):
                break
        if len(independent_indices) < len(A):
            return A[independent_indices], b[independent_indices]
        return A, b
            
    def build_equation_system(self, use_current_betas=True) -> Tuple[np.ndarray, np.ndarray, List]:
        """
        Build the system of linear equations for moles of each component in each tube.
        
        For reactors with stoichiometry aA + bB -> cC + dD and conversion α:
        - Find the limiting reactant (the one with least moles relative to stoichiometry)
        - Express all consumption/production in terms of the limiting reactant consumed
        - For component C: n_out[C] = n_in[C] + (stoich[C]/stoich[lim]) * α * n_in[lim]
        
        Parameters:
            use_current_betas: If True, use current beta values. If False, uses variable beta tubes'
                               solved_beta if available, otherwise their beta value.
        """
        component_ids = list(self.components.keys())
        tube_ids = list(self.tubes.keys())
        
        # Variables: moles of each component in each tube (no extent variables)
        variables = []
        for tid in tube_ids:
            for cid in component_ids:
                variables.append(('tube', tid, cid))
        
        n_vars = len(variables)
        var_index = {v: i for i, v in enumerate(variables)}
        
        equations = []
        
        for proc_id, process in self.processes.items():
            incoming = [t for t in self.tubes.values() if t.to_process == proc_id]
            outgoing = [t for t in self.tubes.values() if t.from_process == proc_id]
            
            if process.process_type == ProcessType.INPUT:
                # Input: outgoing moles = input moles for each component
                for comp_id in component_ids:
                    input_mol = self.input_moles.get(proc_id, {}).get(comp_id, 0.0)
                    if input_mol > 0 or outgoing:
                        coeffs = np.zeros(n_vars)
                        for tube in outgoing:
                            coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                        equations.append((coeffs, input_mol))
                        
            elif process.process_type == ProcessType.REACTOR:
                # Reactor with stoichiometry
                # α = fraction of reference reactant converted
                # n_out[C] = n_in[C] + stoich[C] × α × n_in[ref] / |stoich[ref]|
                # Since stoich[ref] is negative, this correctly handles signs
                conversion = process.conversion
                
                # Find the reference reactant: choose the one with smallest |stoich|
                # E.g., for 3H2 + N2 -> 2NH3, choose N2 (|stoich|=1)
                reference_reactant = None
                reference_stoich = None
                for cid, coeff in process.stoich.items():
                    if coeff < 0:  # This is a reactant
                        if reference_stoich is None or abs(coeff) < abs(reference_stoich):
                            reference_reactant = cid
                            reference_stoich = coeff  # Negative value
                
                for comp_id in component_ids:
                    stoich_coeff = process.stoich.get(comp_id, 0.0)
                    
                    if outgoing and reference_reactant is not None:
                        coeffs = np.zeros(n_vars)
                        
                        # Extent of reaction: ξ = α × n_in[ref] / |stoich[ref]|
                        # Change for component C: Δn[C] = stoich[C] × ξ
                        # So: n_out[C] = n_in[C] + stoich[C] × α × n_in[ref] / |stoich[ref]|
                        #
                        # For ref (stoich=-1): n_out = n_in + (-1) × α × n_in / 1 = n_in × (1-α) ✓
                        # For H2 (stoich=-3): n_out = n_in[H2] + (-3) × α × n_in[N2] / 1 = n_in[H2] - 3α×n_in[N2] ✓
                        # For NH3 (stoich=+2): n_out = n_in[NH3] + 2 × α × n_in[N2] / 1 = n_in[NH3] + 2α×n_in[N2] ✓
                        
                        # n_out[C]
                        for tube in outgoing:
                            coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                        
                        # - n_in[C]
                        for tube in incoming:
                            coeffs[var_index[('tube', tube.id, comp_id)]] = -1.0
                        
                        # - stoich[C] × α / |stoich[ref]| × n_in[ref]
                        if stoich_coeff != 0 and reference_stoich is not None:
                            factor = stoich_coeff * conversion / abs(reference_stoich)
                            for tube in incoming:
                                coeffs[var_index[('tube', tube.id, reference_reactant)]] -= factor
                        
                        equations.append((coeffs, 0.0))
                        
            elif process.process_type == ProcessType.MIXER:
                # Mixer: output = sum of inputs for each component (mass balance)
                for comp_id in component_ids:
                    if outgoing:
                        coeffs = np.zeros(n_vars)
                        for tube in outgoing:
                            coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                        for tube in incoming:
                            coeffs[var_index[('tube', tube.id, comp_id)]] = -1.0
                        equations.append((coeffs, 0.0))
                        
            elif process.process_type == ProcessType.SEPARATOR:
                # Separator: physical split where all components split with the SAME β
                if self.alpha_eff > 0:
                    # Find reference reactant and stoichiometry
                    ref_reactant = None
                    ref_stoich = None
                    all_stoich = {}
                    for proc in self.processes.values():
                        if proc.process_type == ProcessType.REACTOR:
                            for cid, coeff in proc.stoich.items():
                                all_stoich[cid] = coeff
                                if coeff < 0:
                                    if ref_stoich is None or abs(coeff) < abs(ref_stoich):
                                        ref_reactant = cid
                                        ref_stoich = coeff
                    
                    # Get original input values for ratio calculation
                    input_procs = [p for p in self.processes.values() if p.process_type == ProcessType.INPUT]
                    total_input = {}
                    for comp_id in component_ids:
                        total_input[comp_id] = sum(
                            self.input_moles.get(p.id, {}).get(comp_id, 0.0) for p in input_procs
                        )
                    ref_input = total_input.get(ref_reactant, 0.0)
                    
                    # Mass balance for reference component
                    if ref_reactant:
                        coeffs = np.zeros(n_vars)
                        for tube in outgoing:
                            coeffs[var_index[('tube', tube.id, ref_reactant)]] = 1.0
                        for in_tube in incoming:
                            coeffs[var_index[('tube', in_tube.id, ref_reactant)]] = -1.0
                        equations.append((coeffs, 0.0))
                    
                    # For each output tube, constrain component ratios
                    # The key insight: purge[C]/purge[ref] ratio is KNOWN from alpha_eff
                    # And recycle[C]/recycle[ref] must be the SAME ratio (physical splitter)
                    # 
                    # For reactants: ratio = input[C]/input[ref] (stoichiometric feed ratio)
                    # For products: ratio = (stoich[C] * α_eff) / ((1-α_eff) * |stoich[ref]|)
                    #               from purge[prod] = stoich * α_eff * input[ref] / |stoich[ref]|
                    #               and  purge[ref] = input[ref] * (1-α_eff)
                    
                    for tube in outgoing:
                        for comp_id in component_ids:
                            if comp_id != ref_reactant and ref_stoich is not None:
                                stoich_c = all_stoich.get(comp_id, 0.0)
                                
                                if stoich_c < 0:
                                    # Reactant: calculate based on stoichiometry
                                    # out[C] = input[C] - |stoich[C]|/|stoich[ref]| * α_eff * input[ref]
                                    # out[ref] = input[ref] * (1 - α_eff)
                                    # ratio = out[C] / out[ref]
                                    input_c = total_input.get(comp_id, 0.0)
                                    consumed_c = abs(stoich_c) / abs(ref_stoich) * self.alpha_eff * ref_input
                                    out_c = input_c - consumed_c
                                    out_ref = ref_input * (1 - self.alpha_eff)
                                    if out_ref > 0:
                                        ratio = out_c / out_ref
                                    else:
                                        ratio = 0.0
                                elif stoich_c > 0:
                                    # Product: ratio from alpha_eff relationship
                                    # purge[C] = stoich[C] * α_eff * input[ref] / |stoich[ref]|
                                    # purge[ref] = input[ref] * (1 - α_eff)
                                    # ratio = purge[C] / purge[ref] = stoich[C] * α_eff / (|stoich[ref]| * (1-α_eff))
                                    ratio = stoich_c * self.alpha_eff / (abs(ref_stoich) * (1 - self.alpha_eff))
                                else:
                                    # Inert: ratio from input
                                    if ref_input > 0:
                                        ratio = total_input.get(comp_id, 0.0) / ref_input
                                    else:
                                        ratio = 0.0
                                
                                # Constraint: out[C] = ratio * out[ref]
                                coeffs = np.zeros(n_vars)
                                coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                                coeffs[var_index[('tube', tube.id, ref_reactant)]] = -ratio
                                equations.append((coeffs, 0.0))
                else:
                    # No alpha_eff: use beta values to determine split
                    # Key insight: when two components have the SAME beta in a tube,
                    # their output ratio equals their FEED ratio (physical splitter preserves ratios)
                    
                    # Get feed ratios from INPUT processes
                    input_procs = [p for p in self.processes.values() if p.process_type == ProcessType.INPUT]
                    total_feed = {}
                    for comp_id in component_ids:
                        total_feed[comp_id] = sum(
                            self.input_moles.get(inp.id, {}).get(comp_id, 0.0) for inp in input_procs
                        )
                    
                    # First: mass balance for each component (sum of outputs = sum of inputs)
                    for comp_id in component_ids:
                        coeffs = np.zeros(n_vars)
                        for tube in outgoing:
                            coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                        for in_tube in incoming:
                            coeffs[var_index[('tube', in_tube.id, comp_id)]] = -1.0
                        equations.append((coeffs, 0.0))
                    
                    # Second: for each tube, get component betas and add constraints
                    for tube in outgoing:
                        # Get beta for each component
                        comp_betas = {}
                        for comp_id in component_ids:
                            if tube.id in process.beta_components:
                                comp_betas[comp_id] = process.beta_components[tube.id].get(comp_id, 0.0)
                            else:
                                comp_betas[comp_id] = tube.beta
                        
                        # Group components by their beta value
                        beta_groups = {}
                        for comp_id, beta_c in comp_betas.items():
                            if beta_c not in beta_groups:
                                beta_groups[beta_c] = []
                            beta_groups[beta_c].append(comp_id)
                        
                        # For each beta group with beta > 0
                        for beta_c, comps in beta_groups.items():
                            if beta_c > 0 and len(comps) >= 1:
                                # Pick reference component (prefer one with feed > 0)
                                ref_comp = comps[0]
                                for c in comps:
                                    if total_feed.get(c, 0) > 0:
                                        ref_comp = c
                                        break
                                
                                ref_feed = total_feed.get(ref_comp, 0.0)
                                
                                # Add constraint: out[ref] = beta * in[ref]
                                coeffs = np.zeros(n_vars)
                                coeffs[var_index[('tube', tube.id, ref_comp)]] = 1.0
                                for in_tube in incoming:
                                    coeffs[var_index[('tube', in_tube.id, ref_comp)]] = -beta_c
                                equations.append((coeffs, 0.0))
                                
                                # For other components in same beta group:
                                # Use FEED ratio to constrain output ratio
                                # out[C] / out[ref] = feed[C] / feed[ref]
                                # Linear form: out[C] = out[ref] * (feed[C] / feed[ref])
                                for comp_id in comps:
                                    if comp_id != ref_comp:
                                        comp_feed = total_feed.get(comp_id, 0.0)
                                        if ref_feed > 0:
                                            ratio = comp_feed / ref_feed
                                        else:
                                            ratio = 0.0
                                        
                                        # out[C] - ratio * out[ref] = 0
                                        coeffs = np.zeros(n_vars)
                                        coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                                        coeffs[var_index[('tube', tube.id, ref_comp)]] = -ratio
                                        equations.append((coeffs, 0.0))
                        
            elif process.process_type == ProcessType.OUTPUT:
                # Output: just receives flow, no equations needed
                pass
        
        # Add alpha_eff constraint if enabled (> 0)
        if self.alpha_eff > 0:
            # Find input and output processes
            input_procs = [p for p in self.processes.values() if p.process_type == ProcessType.INPUT]
            output_procs = [p for p in self.processes.values() if p.process_type == ProcessType.OUTPUT]
            
            # Find the reference reactant (smallest |stoich| among reactants)
            # and collect stoichiometry
            reference_reactant = None
            reference_stoich = None
            all_stoich = {}
            for proc in self.processes.values():
                if proc.process_type == ProcessType.REACTOR:
                    for cid, coeff in proc.stoich.items():
                        all_stoich[cid] = coeff
                        if coeff < 0:  # Reactant
                            if reference_stoich is None or abs(coeff) < abs(reference_stoich):
                                reference_reactant = cid
                                reference_stoich = coeff
            
            if input_procs and output_procs and reference_reactant:
                # Get input moles for ALL components
                total_input = {}
                for comp_id in component_ids:
                    total_input[comp_id] = sum(
                        self.input_moles.get(inp.id, {}).get(comp_id, 0.0) for inp in input_procs
                    )
                
                ref_total = total_input.get(reference_reactant, 0.0)
                
                if ref_total > 0:
                    # For each component, constrain the output based on alpha_eff
                    for comp_id in component_ids:
                        stoich_c = all_stoich.get(comp_id, 0.0)
                        input_c = total_input.get(comp_id, 0.0)
                        
                        coeffs = np.zeros(n_vars)
                        for out_proc in output_procs:
                            incoming_to_output = [t for t in self.tubes.values() if t.to_process == out_proc.id]
                            for tube in incoming_to_output:
                                coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                        
                        if stoich_c < 0:
                            # Reactant: output = input - consumed
                            # consumed = |stoich[C]| / |stoich[ref]| * α_eff * input[ref]
                            consumed = abs(stoich_c) / abs(reference_stoich) * self.alpha_eff * ref_total
                            rhs = input_c - consumed
                        elif stoich_c > 0:
                            # Product: output = input + produced
                            # produced = ref_consumed * (stoich_product / |stoich_ref|)
                            # ref_consumed = ref_input * α_eff
                            produced = ref_total * self.alpha_eff * stoich_c / abs(reference_stoich)
                            rhs = input_c + produced
                        else:
                            # Inert (stoich = 0): output = input (no conversion)
                            rhs = input_c
                        
                        equations.append((coeffs, rhs))
                    
        if not equations:
            return np.array([]), np.array([]), variables
            
        A = np.array([eq[0] for eq in equations])
        b = np.array([eq[1] for eq in equations])
        
        return A, b, variables
        
    def display_results(self):
        """Display moles of each component using collapsible sections"""
        # Clear existing widgets in results frame
        for widget in self.results_scroll_frame.winfo_children():
            widget.destroy()
        
        # Alpha_eff section
        if self.alpha_eff > 0:
            def alpha_content(parent):
                ttk.Label(parent, text=f"{self.alpha_eff*100:.1f}% of limiting reactant converted").pack(anchor='w')
                ttk.Label(parent, text="Separator splits from α_eff constraint").pack(anchor='w')
            self.create_collapsible_section(
                self.results_scroll_frame, 
                f"α_eff: {self.alpha_eff:.4f} (enabled)",
                alpha_content, 
                "alpha_eff",
                initially_open=False
            )
        else:
            ttk.Label(self.results_scroll_frame, text="α_eff: disabled (using β values)", 
                     font=('Arial', 9)).pack(anchor='w', padx=5, pady=2)
        
        ttk.Separator(self.results_scroll_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=5)
        
        # Reactors section
        reactors = [p for p in self.processes.values() if p.process_type == ProcessType.REACTOR]
        if reactors:
            ttk.Label(self.results_scroll_frame, text="REACTORS", 
                     font=('Arial', 10, 'bold')).pack(anchor='w', padx=5, pady=2)
            
            for proc in reactors:
                conv = proc.conversion
                
                def make_reactor_content(p=proc):
                    def reactor_content(parent):
                        ttk.Label(parent, text=f"Per-pass conversion: {p.conversion:.1%}").pack(anchor='w')
                        
                        if p.stoich:
                            ttk.Label(parent, text="Stoichiometry:", font=('Arial', 8, 'bold')).pack(anchor='w', pady=(5,0))
                            reactants = []
                            products = []
                            for comp_id, coeff in p.stoich.items():
                                comp_name = self.components.get(comp_id, Component(comp_id, comp_id)).name
                                if coeff < 0:
                                    reactants.append(f"{abs(coeff):.0f}{comp_name}")
                                elif coeff > 0:
                                    products.append(f"{coeff:.0f}{comp_name}")
                            
                            if reactants and products:
                                reaction_str = " + ".join(reactants) + " → " + " + ".join(products)
                                ttk.Label(parent, text=reaction_str, font=('Arial', 9)).pack(anchor='w', padx=10)
                    return reactor_content
                
                self.create_collapsible_section(
                    self.results_scroll_frame,
                    f"{proc.name} (α = {conv:.1%})",
                    make_reactor_content(),
                    f"reactor_{proc.id}",
                    initially_open=False
                )
            
            ttk.Separator(self.results_scroll_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=5)
        
        # Tubes section
        if self.tubes:
            ttk.Label(self.results_scroll_frame, text="TUBES", 
                     font=('Arial', 10, 'bold')).pack(anchor='w', padx=5, pady=2)
            
            for tid, tube in self.tubes.items():
                from_name = self.processes[tube.from_process].name
                to_name = self.processes[tube.to_process].name
                total = tube.total_moles() if tube.moles else 0
                display_total = abs(total) if abs(total) < 1e-9 else total
                
                def make_tube_content(t=tube, total_val=display_total):
                    def tube_content(parent):
                        ttk.Label(parent, text=f"β = {t.beta:.3f}").pack(anchor='w')
                        
                        if t.moles:
                            for comp_id, moles in t.moles.items():
                                comp_name = self.components.get(comp_id, Component(comp_id, comp_id)).name
                                display_moles = abs(moles) if abs(moles) < 1e-9 else moles
                                ttk.Label(parent, text=f"{comp_name}: {display_moles:.4f} mol").pack(anchor='w')
                            
                            ttk.Separator(parent, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=2)
                            ttk.Label(parent, text=f"TOTAL: {total_val:.4f} mol", 
                                     font=('Arial', 9, 'bold')).pack(anchor='w')
                    return tube_content
                
                self.create_collapsible_section(
                    self.results_scroll_frame,
                    f"{tid}: {from_name} → {to_name} ({display_total:.2f} mol)",
                    make_tube_content(),
                    f"tube_{tid}",
                    initially_open=False
                )
        
        # Update scroll region
        self.results_scroll_frame.update_idletasks()
        self.results_canvas.configure(scrollregion=self.results_canvas.bbox("all"))
            
    def show_symbolic_equations(self):
        """Show only the symbolic equations used to solve the system"""
        if not self.components:
            messagebox.showinfo("Info", "Add components first", parent=self.root)
            return
            
        # Create new window
        eq_window = tk.Toplevel(self.root)
        eq_window.title("Symbolic Equations")
        eq_window.geometry("700x500")
        eq_window.transient(self.root)
        
        text = tk.Text(eq_window, font=('Consolas', 11))
        scrollbar = ttk.Scrollbar(eq_window, orient='vertical', command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Find reference reactant info
        reference_reactant = None
        reference_stoich = None
        all_stoich = {}
        for proc in self.processes.values():
            if proc.process_type == ProcessType.REACTOR:
                for cid, coeff in proc.stoich.items():
                    all_stoich[cid] = coeff
                    if coeff < 0:
                        if reference_stoich is None or abs(coeff) < abs(reference_stoich):
                            reference_reactant = cid
                            reference_stoich = coeff
        
        eq_num = 1
        component_ids = list(self.components.keys())
        
        for pid, process in self.processes.items():
            incoming = [t for t in self.tubes.values() if t.to_process == pid]
            outgoing = [t for t in self.tubes.values() if t.from_process == pid]
            
            if process.process_type == ProcessType.INPUT:
                for tube in outgoing:
                    for cid in component_ids:
                        cname = self.components.get(cid, Component(cid, cid)).name
                        text.insert(tk.END, f"({eq_num}) n({cname},{tube.id}) = n₀({cname})\n")
                        eq_num += 1
                
            elif process.process_type == ProcessType.MIXER:
                for cid in component_ids:
                    cname = self.components.get(cid, Component(cid, cid)).name
                    out_tubes = " + ".join([f"n({cname},{t.id})" for t in outgoing])
                    in_tubes = " + ".join([f"n({cname},{t.id})" for t in incoming])
                    text.insert(tk.END, f"({eq_num}) {out_tubes} = {in_tubes}\n")
                    eq_num += 1
                
            elif process.process_type == ProcessType.REACTOR:
                # Find local reference
                local_ref = None
                local_ref_stoich = None
                for cid, coeff in process.stoich.items():
                    if coeff < 0:
                        if local_ref_stoich is None or abs(coeff) < abs(local_ref_stoich):
                            local_ref = cid
                            local_ref_stoich = coeff
                
                if local_ref:
                    ref_name = self.components.get(local_ref, Component(local_ref, local_ref)).name
                    
                    for cid in component_ids:
                        cname = self.components.get(cid, Component(cid, cid)).name
                        stoich_c = process.stoich.get(cid, 0.0)
                        
                        out_str = " + ".join([f"n({cname},{t.id})" for t in outgoing])
                        in_str = " + ".join([f"n({cname},{t.id})" for t in incoming])
                        ref_in_str = " + ".join([f"n({ref_name},{t.id})" for t in incoming])
                        
                        if stoich_c < 0:
                            text.insert(tk.END, f"({eq_num}) {out_str} = {in_str} - (|ν_{cname}|/|ν_{ref_name}|)·α·{ref_in_str}\n")
                        elif stoich_c > 0:
                            text.insert(tk.END, f"({eq_num}) {out_str} = {in_str} + (ν_{cname}/|ν_{ref_name}|)·α·{ref_in_str}\n")
                        else:
                            text.insert(tk.END, f"({eq_num}) {out_str} = {in_str}\n")
                        eq_num += 1
                
            elif process.process_type == ProcessType.SEPARATOR:
                if self.alpha_eff > 0 and reference_reactant:
                    ref_name = self.components.get(reference_reactant, Component(reference_reactant, reference_reactant)).name
                    
                    # Mass balance for reference
                    out_str = " + ".join([f"n({ref_name},{t.id})" for t in outgoing])
                    in_str = " + ".join([f"n({ref_name},{t.id})" for t in incoming])
                    text.insert(tk.END, f"({eq_num}) {out_str} = {in_str}\n")
                    eq_num += 1
                    
                    # Ratio constraints
                    for tube in outgoing:
                        for cid in component_ids:
                            if cid != reference_reactant:
                                cname = self.components.get(cid, Component(cid, cid)).name
                                stoich_c = all_stoich.get(cid, 0.0)
                                
                                if stoich_c < 0:
                                    text.insert(tk.END, f"({eq_num}) n({cname},{tube.id}) = (n_out({cname})/n_out({ref_name}))·n({ref_name},{tube.id})\n")
                                elif stoich_c > 0:
                                    text.insert(tk.END, f"({eq_num}) n({cname},{tube.id}) = (ν_{cname}·α_eff/(|ν_{ref_name}|·(1-α_eff)))·n({ref_name},{tube.id})\n")
                                else:
                                    text.insert(tk.END, f"({eq_num}) n({cname},{tube.id}) = (n₀({cname})/n₀({ref_name}))·n({ref_name},{tube.id})\n")
                                eq_num += 1
                else:
                    for cid in component_ids:
                        cname = self.components.get(cid, Component(cid, cid)).name
                        out_str = " + ".join([f"n({cname},{t.id})" for t in outgoing])
                        in_str = " + ".join([f"n({cname},{t.id})" for t in incoming])
                        text.insert(tk.END, f"({eq_num}) {out_str} = {in_str}\n")
                        eq_num += 1
                    
                    for tube in outgoing:
                        for cid in component_ids:
                            cname = self.components.get(cid, Component(cid, cid)).name
                            in_str = " + ".join([f"n({cname},{t.id})" for t in incoming])
                            text.insert(tk.END, f"({eq_num}) n({cname},{tube.id}) = β_{cname}·({in_str})\n")
                            eq_num += 1
        
        text.insert(tk.END, "\n" + "="*60 + "\n")
        text.insert(tk.END, f"Total equations: {eq_num - 1}\n")
            
    def show_equations(self):
        """Show the system of equations in readable, descriptive format"""
        if not self.components:
            messagebox.showinfo("Info", "Add components first", parent=self.root)
            return
            
        # Create new window
        eq_window = tk.Toplevel(self.root)
        eq_window.title("System of Equations")
        eq_window.geometry("900x700")
        eq_window.transient(self.root)
        
        text = tk.Text(eq_window, font=('Consolas', 10))
        scrollbar = ttk.Scrollbar(eq_window, orient='vertical', command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Header
        text.insert(tk.END, "CHEMICAL FLOW NETWORK - EQUATION SYSTEM\n")
        text.insert(tk.END, "="*70 + "\n\n")
        
        # Show settings
        text.insert(tk.END, "SETTINGS:\n")
        text.insert(tk.END, "-"*40 + "\n")
        if self.alpha_eff > 0:
            text.insert(tk.END, f"  α_eff = {self.alpha_eff:.4f} (overall conversion enabled)\n")
        else:
            text.insert(tk.END, "  α_eff = 0 (disabled, using β values for separators)\n")
        text.insert(tk.END, "\n")
        
        # Find reference reactant
        reference_reactant = None
        reference_stoich = None
        all_stoich = {}
        for proc in self.processes.values():
            if proc.process_type == ProcessType.REACTOR:
                for cid, coeff in proc.stoich.items():
                    all_stoich[cid] = coeff
                    if coeff < 0:
                        if reference_stoich is None or abs(coeff) < abs(reference_stoich):
                            reference_reactant = cid
                            reference_stoich = coeff
        
        if reference_reactant:
            ref_name = self.components.get(reference_reactant, Component(reference_reactant, reference_reactant)).name
            text.insert(tk.END, f"  Reference reactant: {ref_name} (stoich = {reference_stoich})\n\n")
        
        # Show input values
        text.insert(tk.END, "INPUT VALUES (known):\n")
        text.insert(tk.END, "-"*40 + "\n")
        input_procs = [p for p in self.processes.values() if p.process_type == ProcessType.INPUT]
        for inp in input_procs:
            if inp.id in self.input_moles:
                for cid, moles in self.input_moles[inp.id].items():
                    cname = self.components.get(cid, Component(cid, cid)).name
                    text.insert(tk.END, f"  n({cname}, from {inp.name}) = {moles:.4f} mol\n")
        text.insert(tk.END, "\n")
        
        eq_num = 1
        
        # Process-by-process equations
        text.insert(tk.END, "EQUATIONS BY PROCESS:\n")
        text.insert(tk.END, "="*70 + "\n\n")
        
        component_ids = list(self.components.keys())
        
        for pid, process in self.processes.items():
            incoming = [t for t in self.tubes.values() if t.to_process == pid]
            outgoing = [t for t in self.tubes.values() if t.from_process == pid]
            
            if process.process_type == ProcessType.INPUT:
                text.insert(tk.END, f"▸ {process.name} (INPUT)\n")
                text.insert(tk.END, "  Known input moles flow to output tube(s)\n")
                for tube in outgoing:
                    to_name = self.processes[tube.to_process].name
                    for cid in component_ids:
                        cname = self.components.get(cid, Component(cid, cid)).name
                        moles = self.input_moles.get(pid, {}).get(cid, 0.0)
                        text.insert(tk.END, f"  ({eq_num}) n({cname}, {tube.id}) = {moles:.4f}\n")
                        eq_num += 1
                text.insert(tk.END, "\n")
                
            elif process.process_type == ProcessType.MIXER:
                text.insert(tk.END, f"▸ {process.name} (MIXER)\n")
                text.insert(tk.END, "  Mass balance: output = sum of inputs\n")
                for cid in component_ids:
                    cname = self.components.get(cid, Component(cid, cid)).name
                    out_tubes = ", ".join([f"n({cname},{t.id})" for t in outgoing])
                    in_tubes = " + ".join([f"n({cname},{t.id})" for t in incoming])
                    text.insert(tk.END, f"  ({eq_num}) {out_tubes} = {in_tubes}\n")
                    eq_num += 1
                text.insert(tk.END, "\n")
                
            elif process.process_type == ProcessType.REACTOR:
                text.insert(tk.END, f"▸ {process.name} (REACTOR, α = {process.conversion:.1%})\n")
                
                # Show reaction
                reactants = [f"{abs(v):.0f}{self.components.get(k, Component(k,k)).name}" 
                            for k, v in process.stoich.items() if v < 0]
                products = [f"{v:.0f}{self.components.get(k, Component(k,k)).name}" 
                           for k, v in process.stoich.items() if v > 0]
                rxn_str = " + ".join(reactants) + " → " + " + ".join(products)
                text.insert(tk.END, f"  Reaction: {rxn_str}\n")
                
                # Find local reference (smallest |stoich|)
                local_ref = None
                local_ref_stoich = None
                for cid, coeff in process.stoich.items():
                    if coeff < 0:
                        if local_ref_stoich is None or abs(coeff) < abs(local_ref_stoich):
                            local_ref = cid
                            local_ref_stoich = coeff
                
                if local_ref:
                    ref_name = self.components.get(local_ref, Component(local_ref, local_ref)).name
                    text.insert(tk.END, f"  Reference: {ref_name} (ν = {local_ref_stoich})\n")
                    text.insert(tk.END, f"  Extent: ξ = α × n_in({ref_name}) / |ν_{ref_name}|\n")
                    text.insert(tk.END, f"  For each component: n_out = n_in + ν × ξ\n\n")
                    
                    for cid in component_ids:
                        cname = self.components.get(cid, Component(cid, cid)).name
                        stoich_c = process.stoich.get(cid, 0.0)
                        
                        out_str = " + ".join([f"n({cname},{t.id})" for t in outgoing])
                        in_str = " + ".join([f"n({cname},{t.id})" for t in incoming])
                        ref_in_str = " + ".join([f"n({ref_name},{t.id})" for t in incoming])
                        
                        if stoich_c != 0:
                            factor = stoich_c * process.conversion / abs(local_ref_stoich)
                            if factor > 0:
                                text.insert(tk.END, f"  ({eq_num}) {out_str} = {in_str} + {factor:.4g} × {ref_in_str}\n")
                            else:
                                text.insert(tk.END, f"  ({eq_num}) {out_str} = {in_str} - {abs(factor):.4g} × {ref_in_str}\n")
                        else:
                            text.insert(tk.END, f"  ({eq_num}) {out_str} = {in_str}  (inert)\n")
                        eq_num += 1
                text.insert(tk.END, "\n")
                
            elif process.process_type == ProcessType.SEPARATOR:
                text.insert(tk.END, f"▸ {process.name} (SEPARATOR)\n")
                
                if self.alpha_eff > 0 and reference_reactant:
                    ref_name = self.components.get(reference_reactant, Component(reference_reactant, reference_reactant)).name
                    text.insert(tk.END, f"  Using α_eff constraint (ratio-based splits)\n")
                    text.insert(tk.END, f"  Mass balance + ratio constraints based on α_eff = {self.alpha_eff:.4f}\n\n")
                    
                    # Mass balance for reference
                    text.insert(tk.END, f"  Mass balance for {ref_name}:\n")
                    out_str = " + ".join([f"n({ref_name},{t.id})" for t in outgoing])
                    in_str = " + ".join([f"n({ref_name},{t.id})" for t in incoming])
                    text.insert(tk.END, f"  ({eq_num}) {out_str} = {in_str}\n")
                    eq_num += 1
                    
                    # Ratio constraints
                    text.insert(tk.END, f"\n  Ratio constraints (each output maintains composition):\n")
                    for tube in outgoing:
                        to_name = self.processes[tube.to_process].name
                        text.insert(tk.END, f"  Tube {tube.id} → {to_name}:\n")
                        for cid in component_ids:
                            if cid != reference_reactant:
                                cname = self.components.get(cid, Component(cid, cid)).name
                                stoich_c = all_stoich.get(cid, 0.0)
                                
                                if stoich_c < 0:
                                    text.insert(tk.END, f"    ({eq_num}) n({cname},{tube.id}) = (out_{cname}/out_{ref_name}) × n({ref_name},{tube.id})\n")
                                elif stoich_c > 0:
                                    ratio = stoich_c * self.alpha_eff / (abs(reference_stoich) * (1 - self.alpha_eff))
                                    text.insert(tk.END, f"    ({eq_num}) n({cname},{tube.id}) = {ratio:.4g} × n({ref_name},{tube.id})  (product ratio)\n")
                                else:
                                    text.insert(tk.END, f"    ({eq_num}) n({cname},{tube.id}) = (in_{cname}/in_{ref_name}) × n({ref_name},{tube.id})  (inert)\n")
                                eq_num += 1
                else:
                    text.insert(tk.END, f"  Using β values for split fractions\n")
                    text.insert(tk.END, f"  For each component: Σ(outputs) = Σ(inputs)\n")
                    text.insert(tk.END, f"  For known β > 0: out = β × in\n\n")
                    
                    for cid in component_ids:
                        cname = self.components.get(cid, Component(cid, cid)).name
                        out_str = " + ".join([f"n({cname},{t.id})" for t in outgoing])
                        in_str = " + ".join([f"n({cname},{t.id})" for t in incoming])
                        text.insert(tk.END, f"  ({eq_num}) {out_str} = {in_str}  (mass balance)\n")
                        eq_num += 1
                    
                    for tube in outgoing:
                        to_name = self.processes[tube.to_process].name
                        text.insert(tk.END, f"\n  Tube {tube.id} → {to_name} (β = {tube.beta:.3f}):\n")
                        for cid in component_ids:
                            cname = self.components.get(cid, Component(cid, cid)).name
                            if tube.id in process.beta_components:
                                beta_c = process.beta_components[tube.id].get(cid, 0.0)
                            else:
                                beta_c = tube.beta
                            if beta_c > 0:
                                in_str = " + ".join([f"n({cname},{t.id})" for t in incoming])
                                text.insert(tk.END, f"    ({eq_num}) n({cname},{tube.id}) = {beta_c:.3f} × ({in_str})\n")
                                eq_num += 1
                text.insert(tk.END, "\n")
                
            elif process.process_type == ProcessType.OUTPUT:
                text.insert(tk.END, f"▸ {process.name} (OUTPUT)\n")
                text.insert(tk.END, "  Terminal node - receives flow, no equations\n\n")
        
        # Alpha_eff output constraints
        if self.alpha_eff > 0:
            text.insert(tk.END, "="*70 + "\n")
            text.insert(tk.END, "OVERALL CONVERSION CONSTRAINTS (α_eff):\n")
            text.insert(tk.END, "-"*40 + "\n")
            ref_name = self.components.get(reference_reactant, Component(reference_reactant, reference_reactant)).name
            text.insert(tk.END, f"Based on α_eff = {self.alpha_eff:.4f} for reference {ref_name}:\n\n")
            
            output_procs = [p for p in self.processes.values() if p.process_type == ProcessType.OUTPUT]
            output_tubes = []
            for out_proc in output_procs:
                output_tubes.extend([t for t in self.tubes.values() if t.to_process == out_proc.id])
            
            input_procs = [p for p in self.processes.values() if p.process_type == ProcessType.INPUT]
            
            for cid in component_ids:
                cname = self.components.get(cid, Component(cid, cid)).name
                stoich_c = all_stoich.get(cid, 0.0)
                
                # Get total input
                total_input = sum(self.input_moles.get(inp.id, {}).get(cid, 0.0) for inp in input_procs)
                ref_total = sum(self.input_moles.get(inp.id, {}).get(reference_reactant, 0.0) for inp in input_procs)
                
                out_str = " + ".join([f"n({cname},{t.id})" for t in output_tubes])
                
                if stoich_c < 0:
                    # Reactant
                    consumed = abs(stoich_c) / abs(reference_stoich) * self.alpha_eff * ref_total
                    expected = total_input - consumed
                    text.insert(tk.END, f"({eq_num}) {out_str} = {total_input:.4f} - {consumed:.4f} = {expected:.4f}\n")
                    text.insert(tk.END, f"     (input - |ν|/|ν_ref| × α_eff × input_ref)\n\n")
                elif stoich_c > 0:
                    # Product
                    produced = ref_total * self.alpha_eff * stoich_c / abs(reference_stoich)
                    expected = total_input + produced
                    text.insert(tk.END, f"({eq_num}) {out_str} = {total_input:.4f} + {produced:.4f} = {expected:.4f}\n")
                    text.insert(tk.END, f"     (input + ν/|ν_ref| × α_eff × input_ref)\n\n")
                else:
                    # Inert
                    text.insert(tk.END, f"({eq_num}) {out_str} = {total_input:.4f}  (inert, no change)\n\n")
                eq_num += 1
            
    def clear_all(self):
        """Clear all processes, tubes and components"""
        if messagebox.askyesno("Confirm", "Clear all components, processes and tubes?", parent=self.root):
            self.save_undo_state()
            self.components.clear()
            self.processes.clear()
            self.tubes.clear()
            self.input_moles.clear()
            self.alpha_eff = 0.0
            self.process_counter = 0
            self.tube_counter = 0
            self.component_counter = 0
            self.redraw()
            self.results_text.delete(1.0, tk.END)
            self.status_var.set("Cleared all")
            
    def save_network(self):
        """Save network to JSON file"""
        from tkinter import filedialog
        
        filename = filedialog.asksaveasfilename(
            parent=self.root,
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            data = {
                "components": {cid: {
                    "id": c.id, "name": c.name
                } for cid, c in self.components.items()},
                "processes": {pid: {
                    "id": p.id, "name": p.name, "type": p.process_type.value,
                    "x": p.x, "y": p.y, "stoich": p.stoich, 
                    "conversion": p.conversion,
                    "beta_components": p.beta_components
                } for pid, p in self.processes.items()},
                "tubes": {tid: {
                    "id": t.id, "from": t.from_process, "to": t.to_process, "beta": t.beta
                } for tid, t in self.tubes.items()},
                "input_moles": self.input_moles,
                "alpha_eff": self.alpha_eff
            }
            
            with open(filename, 'w') as f:
                json.dump(data, f, indent=2)
                
            self.status_var.set(f"Saved to {filename}")
            
    def load_network(self):
        """Load network from JSON file"""
        from tkinter import filedialog
        
        filename = filedialog.askopenfilename(
            parent=self.root,
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            with open(filename, 'r') as f:
                data = json.load(f)
            
            self.components.clear()
            self.processes.clear()
            self.tubes.clear()
            self.input_moles.clear()
            
            # Load components
            for cid, cdata in data.get("components", {}).items():
                self.components[cid] = Component(cdata["id"], cdata["name"])
            
            # Load processes  
            for pid, pdata in data["processes"].items():
                # Support both old 'alpha' and new 'stoich' format
                stoich = pdata.get("stoich", pdata.get("alpha", {}))
                conversion = pdata.get("conversion", 0.0)
                self.processes[pid] = Process(
                    pdata["id"], pdata["name"], ProcessType(pdata["type"]),
                    pdata["x"], pdata["y"], 
                    stoich=stoich,
                    conversion=conversion,
                    beta_components=pdata.get("beta_components", {})
                )
            
            # Load tubes
            for tid, tdata in data["tubes"].items():
                self.tubes[tid] = Tube(
                    tdata["id"], tdata["from"], tdata["to"], 
                    tdata.get("beta", 1.0)
                )
            
            self.input_moles = data.get("input_moles", {})
            self.alpha_eff = data.get("alpha_eff", 0.0)  # Default to disabled
            
            # Update counters based on loaded data
            self.component_counter = len(self.components)
            self.process_counter = len(self.processes)
            self.tube_counter = len(self.tubes)
            
            self.redraw()
            self.status_var.set(f"Loaded from {filename}")
            
    def load_example(self):
        """Load an example network with components"""
        # Don't use clear_all to avoid confirmation dialog
        self.components.clear()
        self.processes.clear()
        self.tubes.clear()
        self.input_moles.clear()
        
        # Define components
        self.components = {
            "A": Component("A", "Reactant A"),
            "B": Component("B", "Reactant B"),
            "C": Component("C", "Product C"),
        }
        
        # Example: Reaction 2A + B -> C with recirculation
        # Stoichiometry: A=-2 (2 mol consumed), B=-1 (1 mol consumed), C=+1 (1 mol produced)
        # Conversion: 0.8 (80% of limiting reactant A converts per pass)
        self.processes = {
            "P1": Process("P1", "Feed", ProcessType.INPUT, 100, 200),
            "P2": Process("P2", "Mixer", ProcessType.MIXER, 250, 200),
            "P3": Process("P3", "Reactor", ProcessType.REACTOR, 400, 200, 
                         stoich={"A": -2, "B": -1, "C": 1},  # 2A + B -> C
                         conversion=0.8),  # 80% conversion per pass
            "P4": Process("P4", "Separator", ProcessType.SEPARATOR, 550, 200,
                         beta_components={
                             "T4": {"A": 0.1, "B": 0.1, "C": 0.95},   # Product stream: mostly C
                             "T5": {"A": 0.9, "B": 0.9, "C": 0.05},   # Recycle: unreacted A, B
                         }),
            "P5": Process("P5", "Product", ProcessType.OUTPUT, 700, 200),
        }
        
        self.tubes = {
            "T1": Tube("T1", "P1", "P2", beta=1.0),
            "T2": Tube("T2", "P2", "P3", beta=1.0),
            "T3": Tube("T3", "P3", "P4", beta=1.0),
            "T4": Tube("T4", "P4", "P5", beta=0.7),   # Default beta (overridden by separator)
            "T5": Tube("T5", "P4", "P2", beta=0.3),   # Recirculation
        }
        
        # Input moles: 100 mol A, 60 mol B, 0 mol C (stoich ratio 2:1, slight excess of B)
        self.input_moles = {"P1": {"A": 100.0, "B": 60.0, "C": 0.0}}
        
        # Overall process efficiency (0 = disabled, only apply if you want to constrain output)
        # If α_eff = 1.0, all of limiting reactant A (100 mol) would be converted
        # So output would have 0 mol A. Set to 0 to disable this constraint.
        self.alpha_eff = 0.0  # Disabled by default
        
        self.process_counter = 5
        self.tube_counter = 5
        self.component_counter = 3  # A, B, C already defined
        
        self.redraw()
        self.status_var.set("Loaded example: 2A + B → C reactor with 80% per-pass conversion (set α_eff to enable overall constraint)")


def main():
    root = tk.Tk()
    app = FlowNetworkGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
