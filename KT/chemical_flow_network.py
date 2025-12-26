"""
Chemical Process Flow Network Calculator

This program models chemical processes connected by tubes/pipes.
- Squares (nodes): Chemical processes, separations, or collections
- Arrows (edges): Tubes that transfer chemicals between processes

Parameters:
- Alpha (α): Transformation factor per component (moles out / moles in for each species)
- Beta (β): Split ratio (fraction going through a specific output tube)
- Alpha-eff: Overall process efficiency (1 - output/input of limiting reactant)
           This applies to the ENTIRE process from system input to output

The program calculates MOLES of each chemical component in ALL TUBES,
solving the system of equations even with recirculation loops.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, field
from enum import Enum


class ProcessType(Enum):
    """Types of chemical processes"""
    INPUT = "input"           # External input source
    OUTPUT = "output"         # External output sink
    REACTOR = "reactor"       # Chemical reaction (uses alpha per component)
    SEPARATOR = "separator"   # Separation process (uses beta for splits, can separate components)
    MIXER = "mixer"           # Collection/mixing point
    SPLITTER = "splitter"     # Flow splitter (uses beta, same composition in all outputs)


@dataclass
class Component:
    """Represents a chemical component/species"""
    id: str
    name: str
    molecular_weight: float = 1.0  # g/mol (optional, for mass calculations)
    
    def __hash__(self):
        return hash(self.id)


@dataclass
class Process:
    """Represents a chemical process (node/square)"""
    id: str
    name: str
    process_type: ProcessType
    # Alpha per component: component_id -> transformation factor
    # For reactors: moles_out = alpha[component] * moles_in
    alpha: Dict[str, float] = field(default_factory=dict)
    # Beta per output tube per component (for separators)
    # tube_id -> component_id -> fraction
    beta_components: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    def __hash__(self):
        return hash(self.id)
    
    def get_alpha(self, component_id: str) -> float:
        """Get alpha for a specific component (default 1.0)"""
        return self.alpha.get(component_id, 1.0)


@dataclass
class Tube:
    """Represents a tube/arrow connecting two processes"""
    id: str
    from_process: str         # Source process ID
    to_process: str           # Destination process ID
    beta: float = 1.0         # Default fraction of source output going through this tube
    # Moles of each component in this tube (solved values)
    moles: Dict[str, float] = field(default_factory=dict)
    
    def __hash__(self):
        return hash(self.id)
    
    def total_moles(self) -> float:
        """Get total moles in tube"""
        return sum(self.moles.values())
    
    def get_moles(self, component_id: str) -> float:
        """Get moles of specific component"""
        return self.moles.get(component_id, 0.0)


@dataclass
class FlowNetwork:
    """Chemical process flow network"""
    components: Dict[str, Component] = field(default_factory=dict)  # All chemical species
    processes: Dict[str, Process] = field(default_factory=dict)
    tubes: Dict[str, Tube] = field(default_factory=dict)
    # Input moles per component: process_id -> component_id -> moles
    input_moles: Dict[str, Dict[str, float]] = field(default_factory=dict)
    # Overall process efficiency (applied to entire system)
    alpha_eff: float = 1.0
    
    def add_component(self, component: Component):
        """Add a chemical component to the system"""
        self.components[component.id] = component
    
    def add_process(self, process: Process):
        """Add a process to the network"""
        self.processes[process.id] = process
        
    def add_tube(self, tube: Tube):
        """Add a tube to the network"""
        self.tubes[tube.id] = tube
    
    def set_alpha_eff(self, alpha_eff: float):
        """Set overall process efficiency"""
        self.alpha_eff = alpha_eff
        
    def set_input_moles(self, process_id: str, component_id: str, moles: float):
        """Set input moles for a specific component at an input process"""
        if process_id not in self.input_moles:
            self.input_moles[process_id] = {}
        self.input_moles[process_id][component_id] = moles
        
    def get_incoming_tubes(self, process_id: str) -> List[Tube]:
        """Get all tubes flowing into a process"""
        return [t for t in self.tubes.values() if t.to_process == process_id]
    
    def get_outgoing_tubes(self, process_id: str) -> List[Tube]:
        """Get all tubes flowing out of a process"""
        return [t for t in self.tubes.values() if t.from_process == process_id]
    
    def get_all_components(self) -> Set[str]:
        """Get all component IDs in the system"""
        return set(self.components.keys())
    
    def build_equation_system(self) -> Tuple[np.ndarray, np.ndarray, List[Tuple[str, str]]]:
        """
        Build the system of linear equations for mass balance.
        Solves for moles of each component in each tube.
        
        Returns:
            A: Coefficient matrix
            b: Constants vector
            variables: List of (tube_id, component_id) pairs (variable order)
        """
        component_ids = list(self.components.keys())
        tube_ids = list(self.tubes.keys())
        
        # Variables: moles of each component in each tube
        # Order: [(T1, C1), (T1, C2), ..., (T2, C1), (T2, C2), ...]
        variables = [(tid, cid) for tid in tube_ids for cid in component_ids]
        n_vars = len(variables)
        var_index = {v: i for i, v in enumerate(variables)}
        
        equations = []  # List of (coefficients, constant)
        
        for proc_id, process in self.processes.items():
            incoming = self.get_incoming_tubes(proc_id)
            outgoing = self.get_outgoing_tubes(proc_id)
            
            # For each component, create mass balance equation
            for comp_id in component_ids:
                
                if process.process_type == ProcessType.INPUT:
                    # Input node: moles out = input moles
                    input_moles = self.input_moles.get(proc_id, {}).get(comp_id, 0.0)
                    if input_moles > 0 or outgoing:
                        coeffs = np.zeros(n_vars)
                        for tube in outgoing:
                            coeffs[var_index[(tube.id, comp_id)]] = 1.0
                        equations.append((coeffs, input_moles))
                        
                elif process.process_type == ProcessType.OUTPUT:
                    # Output node: apply alpha_eff to represent overall efficiency
                    # The output receives incoming * alpha_eff (system-wide)
                    pass  # Handled implicitly; outputs are sinks
                    
                elif process.process_type == ProcessType.REACTOR:
                    # Reactor: output_moles[c] = alpha[c] * input_moles[c]
                    if outgoing:
                        alpha_c = process.get_alpha(comp_id)
                        coeffs = np.zeros(n_vars)
                        for tube in outgoing:
                            coeffs[var_index[(tube.id, comp_id)]] = 1.0
                        for tube in incoming:
                            coeffs[var_index[(tube.id, comp_id)]] = -alpha_c
                        equations.append((coeffs, 0.0))
                        
                elif process.process_type == ProcessType.MIXER:
                    # Mixer: output = sum of all inputs (conservation per component)
                    if outgoing:
                        coeffs = np.zeros(n_vars)
                        for tube in outgoing:
                            coeffs[var_index[(tube.id, comp_id)]] = 1.0
                        for tube in incoming:
                            coeffs[var_index[(tube.id, comp_id)]] = -1.0
                        equations.append((coeffs, 0.0))
                        
                elif process.process_type == ProcessType.SPLITTER:
                    # Splitter: each output gets beta fraction, SAME composition
                    for tube in outgoing:
                        coeffs = np.zeros(n_vars)
                        coeffs[var_index[(tube.id, comp_id)]] = 1.0
                        for in_tube in incoming:
                            coeffs[var_index[(in_tube.id, comp_id)]] = -tube.beta
                        equations.append((coeffs, 0.0))
                        
                elif process.process_type == ProcessType.SEPARATOR:
                    # Separator: can have different beta per component per output
                    for tube in outgoing:
                        coeffs = np.zeros(n_vars)
                        coeffs[var_index[(tube.id, comp_id)]] = 1.0
                        # Get component-specific beta or default to tube beta
                        beta_c = tube.beta
                        if tube.id in process.beta_components:
                            beta_c = process.beta_components[tube.id].get(comp_id, tube.beta)
                        for in_tube in incoming:
                            coeffs[var_index[(in_tube.id, comp_id)]] = -beta_c
                        equations.append((coeffs, 0.0))
        
        # Build matrix
        if not equations:
            return np.array([]), np.array([]), variables
            
        A = np.array([eq[0] for eq in equations])
        b = np.array([eq[1] for eq in equations])
        
        return A, b, variables
    
    def solve(self) -> Dict[str, Dict[str, float]]:
        """
        Solve the flow network and return moles for all components in all tubes.
        
        Returns:
            Dictionary: tube_id -> {component_id: moles}
        """
        A, b, variables = self.build_equation_system()
        
        if len(A) == 0:
            print("No equations to solve!")
            return {}
        
        n_equations, n_variables = A.shape
        
        print(f"\nSystem has {n_equations} equations and {n_variables} variables")
        print(f"  ({len(self.tubes)} tubes × {len(self.components)} components)")
        
        # Remove redundant equations
        A_reduced, b_reduced = self._remove_redundant_equations(A, b)
        n_eq_reduced = len(A_reduced)
        
        if n_eq_reduced < n_variables:
            print(f"Warning: Underdetermined system ({n_eq_reduced} independent equations)")
            solution, residuals, rank, s = np.linalg.lstsq(A_reduced, b_reduced, rcond=None)
        elif n_eq_reduced > n_variables:
            print(f"Overdetermined system. Using least squares.")
            solution, residuals, rank, s = np.linalg.lstsq(A_reduced, b_reduced, rcond=None)
        else:
            try:
                solution = np.linalg.solve(A_reduced, b_reduced)
            except np.linalg.LinAlgError:
                print("Singular matrix. Using least squares.")
                solution, residuals, rank, s = np.linalg.lstsq(A_reduced, b_reduced, rcond=None)
        
        # Apply overall alpha_eff to outputs
        # (scale all flows by alpha_eff for the limiting reactant efficiency)
        
        # Store results in tubes
        results = {tid: {} for tid in self.tubes.keys()}
        for i, (tube_id, comp_id) in enumerate(variables):
            moles = solution[i] * self.alpha_eff  # Apply overall efficiency
            self.tubes[tube_id].moles[comp_id] = moles
            results[tube_id][comp_id] = moles
            
        return results
    
    def _remove_redundant_equations(self, A: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Remove linearly dependent equations"""
        if len(A) == 0:
            return A, b
            
        # Use QR decomposition to find independent rows
        Q, R = np.linalg.qr(A.T)
        independent = np.abs(np.diag(R)) > 1e-10
        
        # Keep only independent rows
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
    
    def print_results(self):
        """Print moles of each component in all tubes"""
        print("\n" + "="*70)
        print("TUBE COMPOSITIONS (MOLES)")
        print("="*70)
        print(f"\nOverall Process Efficiency (α_eff): {self.alpha_eff:.4f}")
        
        # Header
        comp_ids = list(self.components.keys())
        header = f"{'Tube':<10} {'From → To':<25}"
        for cid in comp_ids:
            header += f" {cid:>10}"
        header += f" {'Total':>10}"
        print("\n" + header)
        print("-" * len(header))
        
        for tube in self.tubes.values():
            from_name = self.processes[tube.from_process].name[:10]
            to_name = self.processes[tube.to_process].name[:10]
            
            row = f"{tube.id:<10} {from_name} → {to_name:<12}"
            for cid in comp_ids:
                moles = tube.moles.get(cid, 0.0)
                row += f" {moles:>10.4f}"
            row += f" {tube.total_moles():>10.4f}"
            print(row)
            
        print("\n" + "-"*70)
        print("COMPONENT SUMMARY")
        print("-"*70)
        
        for comp_id, comp in self.components.items():
            print(f"\n{comp.name} ({comp_id}):")
            for tube in self.tubes.values():
                moles = tube.moles.get(comp_id, 0.0)
                if abs(moles) > 1e-10:
                    from_name = self.processes[tube.from_process].name
                    to_name = self.processes[tube.to_process].name
                    print(f"  {tube.id}: {from_name} → {to_name}: {moles:.4f} mol")
    
    def print_equations(self):
        """Print the system of equations for debugging"""
        A, b, variables = self.build_equation_system()
        
        print("\n" + "="*70)
        print("SYSTEM OF EQUATIONS")
        print("="*70)
        print(f"\nVariables ({len(variables)}):")
        for tid in self.tubes.keys():
            tube_vars = [f"{tid}[{cid}]" for cid in self.components.keys()]
            print(f"  {', '.join(tube_vars)}")
        
        print(f"\nEquations ({len(A)}):")
        
        for i, row in enumerate(A):
            terms = []
            for j, coeff in enumerate(row):
                if abs(coeff) > 1e-10:
                    tid, cid = variables[j]
                    var_name = f"{tid}[{cid}]"
                    if coeff == 1.0:
                        terms.append(f"+{var_name}")
                    elif coeff == -1.0:
                        terms.append(f"-{var_name}")
                    elif coeff > 0:
                        terms.append(f"+{coeff:.3f}*{var_name}")
                    else:
                        terms.append(f"{coeff:.3f}*{var_name}")
            equation = " ".join(terms) + f" = {b[i]:.3f}"
            print(f"  Eq{i+1}: {equation}")


def create_example_network() -> FlowNetwork:
    """
    Create an example network with recirculation.
    
    Example: Reactor converting A -> B with recirculation
    
    INPUT(A,B) -> MIXER -> REACTOR -> SEPARATOR -> OUTPUT
                    ^                     |
                    |_____________________|
                    (recirculation of unreacted A)
    """
    network = FlowNetwork()
    
    # Define chemical components
    network.add_component(Component("A", "Reactant A"))
    network.add_component(Component("B", "Product B"))
    
    # Add processes
    network.add_process(Process("P1", "Feed", ProcessType.INPUT))
    network.add_process(Process("P2", "Mixer", ProcessType.MIXER))
    network.add_process(Process(
        "P3", "Reactor", ProcessType.REACTOR,
        alpha={"A": 0.3, "B": 2.5}  # A consumed, B produced (relative to input)
    ))
    network.add_process(Process(
        "P4", "Separator", ProcessType.SEPARATOR,
        beta_components={
            "T4": {"A": 0.1, "B": 0.95},   # To product: 10% A, 95% B
            "T5": {"A": 0.9, "B": 0.05},   # Recirculate: 90% A, 5% B
        }
    ))
    network.add_process(Process("P5", "Product", ProcessType.OUTPUT))
    
    # Add tubes
    network.add_tube(Tube("T1", "P1", "P2", beta=1.0))      # Feed to mixer
    network.add_tube(Tube("T2", "P2", "P3", beta=1.0))      # Mixer to reactor
    network.add_tube(Tube("T3", "P3", "P4", beta=1.0))      # Reactor to separator
    network.add_tube(Tube("T4", "P4", "P5", beta=0.7))      # Separator to output
    network.add_tube(Tube("T5", "P4", "P2", beta=0.3))      # Recirculation
    
    # Set input moles (feed composition)
    network.set_input_moles("P1", "A", 100.0)  # 100 mol/time of A
    network.set_input_moles("P1", "B", 0.0)    # 0 mol/time of B (pure feed)
    
    # Set overall process efficiency
    network.set_alpha_eff(0.95)  # 95% overall efficiency
    
    return network


def create_custom_network() -> FlowNetwork:
    """Interactive network creation"""
    network = FlowNetwork()
    
    print("\n" + "="*70)
    print("CHEMICAL FLOW NETWORK BUILDER")
    print("="*70)
    
    # Add components first
    print("\n--- DEFINE CHEMICAL COMPONENTS ---")
    print("Enter component ID and name. Type 'done' when finished.")
    
    while True:
        comp_id = input("\nComponent ID (or 'done'): ").strip()
        if comp_id.lower() == 'done':
            break
        comp_name = input("Component name: ").strip()
        mw = float(input("Molecular weight (g/mol, default 1.0): ") or 1.0)
        network.add_component(Component(comp_id, comp_name, mw))
        print(f"Added component: {comp_name} ({comp_id})")
    
    if not network.components:
        print("No components defined! Adding default component 'X'")
        network.add_component(Component("X", "Default", 1.0))
    
    # Set overall efficiency
    alpha_eff = float(input("\nOverall process efficiency (α_eff, 0-1, default 1.0): ") or 1.0)
    network.set_alpha_eff(alpha_eff)
    
    # Add processes
    print("\n--- ADD PROCESSES ---")
    print("Process types: input, output, reactor, separator, mixer, splitter")
    print("Enter 'done' when finished adding processes")
    
    while True:
        proc_id = input("\nProcess ID (or 'done'): ").strip()
        if proc_id.lower() == 'done':
            break
            
        name = input("Process name: ").strip()
        ptype = input("Process type: ").strip().lower()
        
        try:
            process_type = ProcessType(ptype)
        except ValueError:
            print(f"Invalid type. Using 'mixer'")
            process_type = ProcessType.MIXER
        
        alpha = {}
        if process_type == ProcessType.REACTOR:
            print("Enter alpha (transformation factor) for each component:")
            for comp_id in network.components:
                a = float(input(f"  Alpha for {comp_id} (default 1.0): ") or 1.0)
                alpha[comp_id] = a
        
        network.add_process(Process(proc_id, name, process_type, alpha=alpha))
        print(f"Added process: {name}")
    
    # Add tubes
    print("\n--- ADD TUBES ---")
    print("Enter 'done' when finished adding tubes")
    
    while True:
        tube_id = input("\nTube ID (or 'done'): ").strip()
        if tube_id.lower() == 'done':
            break
            
        from_proc = input("From process ID: ").strip()
        to_proc = input("To process ID: ").strip()
        beta = float(input("Beta (default split fraction, default 1.0): ") or 1.0)
        
        network.add_tube(Tube(tube_id, from_proc, to_proc, beta))
        print(f"Added tube: {from_proc} -> {to_proc}")
        
        # For separators, ask about component-specific betas
        from_process = network.processes.get(from_proc)
        if from_process and from_process.process_type == ProcessType.SEPARATOR:
            use_comp_beta = input("Set component-specific betas for this tube? (y/n): ").lower()
            if use_comp_beta == 'y':
                if tube_id not in from_process.beta_components:
                    from_process.beta_components[tube_id] = {}
                for comp_id in network.components:
                    b = float(input(f"  Beta for {comp_id} (default {beta}): ") or beta)
                    from_process.beta_components[tube_id][comp_id] = b
    
    # Set input moles
    print("\n--- SET INPUT MOLES ---")
    for proc_id, proc in network.processes.items():
        if proc.process_type == ProcessType.INPUT:
            print(f"\nInput process: {proc.name}")
            for comp_id in network.components:
                moles = float(input(f"  Moles of {comp_id}: ") or 0.0)
                if moles > 0:
                    network.set_input_moles(proc_id, comp_id, moles)
    
    return network


def main():
    print("\n" + "="*70)
    print("CHEMICAL PROCESS FLOW NETWORK CALCULATOR")
    print("Calculates moles of each component in all tubes")
    print("="*70)
    
    print("\nOptions:")
    print("1. Run example network (A -> B reactor with recirculation)")
    print("2. Create custom network")
    
    choice = input("\nSelect option (1 or 2): ").strip()
    
    if choice == "2":
        network = create_custom_network()
    else:
        network = create_example_network()
        print("\nUsing example network: Reactor A->B with recirculation")
    
    # Print the equations
    network.print_equations()
    
    # Solve
    print("\nSolving...")
    results = network.solve()
    
    # Print results
    network.print_results()
    
    return network


if __name__ == "__main__":
    main()
