"""
Chemical Flow Network - Network Model

Contains the data structures that define a chemical flow network:
- ProcessType: Types of process nodes (INPUT, OUTPUT, REACTOR, MIXER, SEPARATOR)
- Component: Chemical species/compounds
- Process: Process nodes in the network
- Tube: Connections between processes
"""

from typing import Dict
from dataclasses import dataclass, field
from enum import Enum


class ProcessType(Enum):
    """Types of chemical processes in the network"""
    INPUT = "input"           # External input source - feeds material into the system
    OUTPUT = "output"         # External output sink - removes material from the system
    REACTOR = "reactor"       # Chemical reaction - transforms reactants to products
    MIXER = "mixer"           # Combines multiple input streams into one output
    SEPARATOR = "separator"   # Splits one input stream into multiple outputs


@dataclass
class Component:
    """
    Chemical component/species in the network.
    
    Attributes:
        id: Unique identifier (e.g., "H2", "N2", "NH3")
        name: Display name (e.g., "Hydrogen", "Nitrogen", "Ammonia")
    """
    id: str
    name: str


@dataclass
class Process:
    """
    Process node in the chemical flow network.
    
    Attributes:
        id: Unique identifier (e.g., "P1", "P2")
        name: Display name (e.g., "Feed", "Reactor", "Separator")
        process_type: Type of process (INPUT, OUTPUT, REACTOR, MIXER, SEPARATOR)
        stoich: Stoichiometric coefficients for reactors
                Negative = reactant (consumed), Positive = product (produced)
                Example: {"N2": -1, "H2": -3, "NH3": 2} for N2 + 3H2 -> 2NH3
        conversion: Fraction of limiting reactant converted per pass (0-1)
                   Only used for REACTOR type
        beta_components: Split fractions for separators
                        Dict of tube_id -> {component_id -> fraction}
                        Example: {"T4": {"NH3": 0.95, "N2": 0.1}}
    """
    id: str
    name: str
    process_type: ProcessType
    stoich: Dict[str, float] = field(default_factory=dict)
    conversion: float = 0.0
    beta_components: Dict[str, Dict[str, float]] = field(default_factory=dict)


@dataclass
class Tube:
    """
    Tube/connection between processes in the network.
    
    Attributes:
        id: Unique identifier (e.g., "T1", "T2")
        from_process: ID of the source process
        to_process: ID of the destination process
        beta: Default split fraction (used when alpha_eff = 0)
        moles: Calculated moles of each component flowing through
               Dict of component_id -> moles
    """
    id: str
    from_process: str
    to_process: str
    beta: float = 1.0
    moles: Dict[str, float] = field(default_factory=dict)
    
    def total_moles(self) -> float:
        """Calculate total moles flowing through this tube"""
        return sum(self.moles.values())
