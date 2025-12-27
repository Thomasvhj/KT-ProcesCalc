"""
Chemical Flow Network - Calculation Module

Contains all mathematical calculations for solving the flow network:
- Building the system of linear equations
- Solving for moles of each component in each tube
- Handling reactors, mixers, separators with alpha/beta parameters
"""

import numpy as np
from typing import Dict, List, Tuple, Optional

# Import network data structures
from chemical_flow_network import ProcessType, Component, Process, Tube


def remove_redundant_equations(A: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Remove linearly dependent equations from the system"""
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


def build_equation_system(
    components: Dict[str, Component],
    processes: Dict[str, Process],
    tubes: Dict[str, Tube],
    input_moles: Dict[str, Dict[str, float]],
    alpha_eff: float
) -> Tuple[np.ndarray, np.ndarray, List]:
    """
    Build the system of linear equations for moles of each component in each tube.
    
    For reactors with stoichiometry aA + bB -> cC + dD and conversion α:
    - Find the reference reactant (the one with smallest |stoichiometry|)
    - Express all consumption/production in terms of the reference reactant
    - For component C: n_out[C] = n_in[C] + (stoich[C]/|stoich[ref]|) * α * n_in[ref]
    
    Parameters:
        components: Dictionary of component_id -> Component
        processes: Dictionary of process_id -> Process
        tubes: Dictionary of tube_id -> Tube
        input_moles: Dictionary of process_id -> {component_id -> moles}
        alpha_eff: Overall conversion efficiency (0 = disabled, use betas)
        
    Returns:
        A: Coefficient matrix
        b: Constants vector
        variables: List of ('tube', tube_id, component_id) tuples
    """
    component_ids = list(components.keys())
    tube_ids = list(tubes.keys())
    
    # Variables: moles of each component in each tube
    variables = []
    for tid in tube_ids:
        for cid in component_ids:
            variables.append(('tube', tid, cid))
    
    n_vars = len(variables)
    var_index = {v: i for i, v in enumerate(variables)}
    
    equations = []
    
    for proc_id, process in processes.items():
        incoming = [t for t in tubes.values() if t.to_process == proc_id]
        outgoing = [t for t in tubes.values() if t.from_process == proc_id]
        
        if process.process_type == ProcessType.INPUT:
            # Input: outgoing moles = input moles for each component
            for comp_id in component_ids:
                input_mol = input_moles.get(proc_id, {}).get(comp_id, 0.0)
                if input_mol > 0 or outgoing:
                    coeffs = np.zeros(n_vars)
                    for tube in outgoing:
                        coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                    equations.append((coeffs, input_mol))
                    
        elif process.process_type == ProcessType.REACTOR:
            # Reactor with stoichiometry
            conversion = process.conversion
            
            # Find the reference reactant: choose the one with smallest |stoich|
            reference_reactant = None
            reference_stoich = None
            for cid, coeff in process.stoich.items():
                if coeff < 0:  # This is a reactant
                    if reference_stoich is None or abs(coeff) < abs(reference_stoich):
                        reference_reactant = cid
                        reference_stoich = coeff
            
            for comp_id in component_ids:
                stoich_coeff = process.stoich.get(comp_id, 0.0)
                
                if outgoing and reference_reactant is not None:
                    coeffs = np.zeros(n_vars)
                    
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
            # Mixer: output = sum of inputs for each component
            for comp_id in component_ids:
                if outgoing:
                    coeffs = np.zeros(n_vars)
                    for tube in outgoing:
                        coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                    for tube in incoming:
                        coeffs[var_index[('tube', tube.id, comp_id)]] = -1.0
                    equations.append((coeffs, 0.0))
                    
        elif process.process_type == ProcessType.SEPARATOR:
            if alpha_eff > 0:
                # Find reference reactant and stoichiometry
                ref_reactant = None
                ref_stoich = None
                all_stoich = {}
                for proc in processes.values():
                    if proc.process_type == ProcessType.REACTOR:
                        for cid, coeff in proc.stoich.items():
                            all_stoich[cid] = coeff
                            if coeff < 0:
                                if ref_stoich is None or abs(coeff) < abs(ref_stoich):
                                    ref_reactant = cid
                                    ref_stoich = coeff
                
                # Get original input values for ratio calculation
                input_procs = [p for p in processes.values() if p.process_type == ProcessType.INPUT]
                total_input = {}
                for comp_id in component_ids:
                    total_input[comp_id] = sum(
                        input_moles.get(p.id, {}).get(comp_id, 0.0) for p in input_procs
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
                for tube in outgoing:
                    for comp_id in component_ids:
                        if comp_id != ref_reactant and ref_stoich is not None:
                            stoich_c = all_stoich.get(comp_id, 0.0)
                            
                            if stoich_c < 0:
                                # Reactant
                                input_c = total_input.get(comp_id, 0.0)
                                consumed_c = abs(stoich_c) / abs(ref_stoich) * alpha_eff * ref_input
                                out_c = input_c - consumed_c
                                out_ref = ref_input * (1 - alpha_eff)
                                if out_ref > 0:
                                    ratio = out_c / out_ref
                                else:
                                    ratio = 0.0
                            elif stoich_c > 0:
                                # Product
                                ratio = stoich_c * alpha_eff / (abs(ref_stoich) * (1 - alpha_eff))
                            else:
                                # Inert
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
                # No alpha_eff: use beta values
                input_procs = [p for p in processes.values() if p.process_type == ProcessType.INPUT]
                total_feed = {}
                for comp_id in component_ids:
                    total_feed[comp_id] = sum(
                        input_moles.get(inp.id, {}).get(comp_id, 0.0) for inp in input_procs
                    )
                
                # Mass balance for each component
                for comp_id in component_ids:
                    coeffs = np.zeros(n_vars)
                    for tube in outgoing:
                        coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                    for in_tube in incoming:
                        coeffs[var_index[('tube', in_tube.id, comp_id)]] = -1.0
                    equations.append((coeffs, 0.0))
                
                # Beta constraints for each tube
                for tube in outgoing:
                    comp_betas = {}
                    for comp_id in component_ids:
                        if tube.id in process.beta_components:
                            comp_betas[comp_id] = process.beta_components[tube.id].get(comp_id, 0.0)
                        else:
                            comp_betas[comp_id] = tube.beta
                    
                    # Group components by beta value
                    beta_groups = {}
                    for comp_id, beta_c in comp_betas.items():
                        if beta_c not in beta_groups:
                            beta_groups[beta_c] = []
                        beta_groups[beta_c].append(comp_id)
                    
                    for beta_c, comps in beta_groups.items():
                        if beta_c > 0 and len(comps) >= 1:
                            # Pick reference component
                            ref_comp = comps[0]
                            for c in comps:
                                if total_feed.get(c, 0) > 0:
                                    ref_comp = c
                                    break
                            
                            ref_feed = total_feed.get(ref_comp, 0.0)
                            
                            # out[ref] = beta * in[ref]
                            coeffs = np.zeros(n_vars)
                            coeffs[var_index[('tube', tube.id, ref_comp)]] = 1.0
                            for in_tube in incoming:
                                coeffs[var_index[('tube', in_tube.id, ref_comp)]] = -beta_c
                            equations.append((coeffs, 0.0))
                            
                            # Ratio constraints for other components in same beta group
                            for comp_id in comps:
                                if comp_id != ref_comp:
                                    comp_feed = total_feed.get(comp_id, 0.0)
                                    if ref_feed > 0:
                                        ratio = comp_feed / ref_feed
                                    else:
                                        ratio = 0.0
                                    
                                    coeffs = np.zeros(n_vars)
                                    coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                                    coeffs[var_index[('tube', tube.id, ref_comp)]] = -ratio
                                    equations.append((coeffs, 0.0))
                    
        elif process.process_type == ProcessType.OUTPUT:
            # Output: just receives flow, no equations needed
            pass
    
    # Add alpha_eff constraint if enabled
    if alpha_eff > 0:
        input_procs = [p for p in processes.values() if p.process_type == ProcessType.INPUT]
        output_procs = [p for p in processes.values() if p.process_type == ProcessType.OUTPUT]
        
        # Find the reference reactant
        reference_reactant = None
        reference_stoich = None
        all_stoich = {}
        for proc in processes.values():
            if proc.process_type == ProcessType.REACTOR:
                for cid, coeff in proc.stoich.items():
                    all_stoich[cid] = coeff
                    if coeff < 0:
                        if reference_stoich is None or abs(coeff) < abs(reference_stoich):
                            reference_reactant = cid
                            reference_stoich = coeff
        
        if input_procs and output_procs and reference_reactant:
            # Get input moles for ALL components
            total_input = {}
            for comp_id in component_ids:
                total_input[comp_id] = sum(
                    input_moles.get(inp.id, {}).get(comp_id, 0.0) for inp in input_procs
                )
            
            ref_total = total_input.get(reference_reactant, 0.0)
            
            if ref_total > 0:
                for comp_id in component_ids:
                    stoich_c = all_stoich.get(comp_id, 0.0)
                    input_c = total_input.get(comp_id, 0.0)
                    
                    coeffs = np.zeros(n_vars)
                    for out_proc in output_procs:
                        incoming_to_output = [t for t in tubes.values() if t.to_process == out_proc.id]
                        for tube in incoming_to_output:
                            coeffs[var_index[('tube', tube.id, comp_id)]] = 1.0
                    
                    if stoich_c < 0:
                        consumed = abs(stoich_c) / abs(reference_stoich) * alpha_eff * ref_total
                        rhs = input_c - consumed
                    elif stoich_c > 0:
                        produced = ref_total * alpha_eff * stoich_c / abs(reference_stoich)
                        rhs = input_c + produced
                    else:
                        rhs = input_c
                    
                    equations.append((coeffs, rhs))
                
    if not equations:
        return np.array([]), np.array([]), variables
        
    A = np.array([eq[0] for eq in equations])
    b = np.array([eq[1] for eq in equations])
    
    return A, b, variables


def solve_network(
    components: Dict[str, Component],
    processes: Dict[str, Process],
    tubes: Dict[str, Tube],
    input_moles: Dict[str, Dict[str, float]],
    alpha_eff: float
) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    """
    Solve the flow network - calculates moles of each component in each tube.
    
    Parameters:
        components: Dictionary of component_id -> Component
        processes: Dictionary of process_id -> Process
        tubes: Dictionary of tube_id -> Tube
        input_moles: Dictionary of process_id -> {component_id -> moles}
        alpha_eff: Overall conversion efficiency (0 = disabled)
        
    Returns:
        results: Dictionary of tube_id -> {component_id -> moles}
        warnings: List of warning messages
    """
    warnings = []
    
    # Build equation system
    A, b, variables = build_equation_system(
        components, processes, tubes, input_moles, alpha_eff
    )
    
    if len(A) == 0:
        return {}, ["No equations to solve!"]
        
    n_equations, n_variables = A.shape
    
    # Remove redundant equations
    A_reduced, b_reduced = remove_redundant_equations(A, b)
    
    # Solve
    if len(A_reduced) == n_variables:
        try:
            solution = np.linalg.solve(A_reduced, b_reduced)
        except np.linalg.LinAlgError:
            solution, _, _, _ = np.linalg.lstsq(A_reduced, b_reduced, rcond=None)
    else:
        solution, _, _, _ = np.linalg.lstsq(A_reduced, b_reduced, rcond=None)
    
    # Check for negative mole values
    negative_vars = []
    for i, var in enumerate(variables):
        if var[0] == 'tube' and solution[i] < -1e-6:
            negative_vars.append((var[1], var[2], solution[i]))
    
    if negative_vars:
        msg = "Some intermediate flows have negative values:\n"
        for tid, cid, val in negative_vars[:5]:
            msg += f"  {tid} {cid}: {val:.4f}\n"
        if len(negative_vars) > 5:
            msg += f"  ... and {len(negative_vars) - 5} more\n"
        warnings.append(msg)
    
    # Build results dictionary
    results = {}
    for i, var in enumerate(variables):
        var_type = var[0]
        if var_type == 'tube':
            _, tube_id, comp_id = var
            if tube_id not in results:
                results[tube_id] = {}
            results[tube_id][comp_id] = solution[i]
    
    return results, warnings


def validate_separator_betas(
    processes: Dict[str, Process],
    tubes: Dict[str, Tube],
    components: Dict[str, Component]
) -> List[str]:
    """
    Validate that separator betas sum to 1 for each component.
    
    Returns:
        List of warning messages (empty if valid)
    """
    warnings = []
    
    for proc in processes.values():
        if proc.process_type == ProcessType.SEPARATOR:
            outgoing = [t for t in tubes.values() if t.from_process == proc.id]
            if len(outgoing) > 1:
                for comp_id in components.keys():
                    beta_sum = 0.0
                    for tube in outgoing:
                        if tube.id in proc.beta_components:
                            beta_sum += proc.beta_components[tube.id].get(comp_id, tube.beta)
                        else:
                            beta_sum += tube.beta
                    if abs(beta_sum - 1.0) > 0.01:
                        comp_name = components.get(comp_id, Component(comp_id, comp_id)).name
                        warnings.append(f"{proc.name}: β sum for {comp_name} = {beta_sum:.2f} (should be 1.0)")
    
    return warnings


def find_reference_reactant(processes: Dict[str, Process]) -> Tuple[Optional[str], Optional[float], Dict[str, float]]:
    """
    Find the reference reactant (smallest |stoichiometry| among reactants).
    
    Returns:
        reference_id: Component ID of reference reactant (or None)
        reference_stoich: Stoichiometric coefficient (negative, or None)
        all_stoich: Dictionary of all stoichiometric coefficients
    """
    reference_reactant = None
    reference_stoich = None
    all_stoich = {}
    
    for proc in processes.values():
        if proc.process_type == ProcessType.REACTOR:
            for cid, coeff in proc.stoich.items():
                all_stoich[cid] = coeff
                if coeff < 0:
                    if reference_stoich is None or abs(coeff) < abs(reference_stoich):
                        reference_reactant = cid
                        reference_stoich = coeff
    
    return reference_reactant, reference_stoich, all_stoich


def calculate_stoichiometric_input_ratios(
    processes: Dict[str, Process],
    components: Dict[str, Component]
) -> Dict[str, float]:
    """
    Calculate stoichiometric input ratios based on reactor stoichiometry.
    
    For reactants (negative stoichiometry): ratio = |stoich| / |reference_stoich|
    For products (positive stoichiometry): ratio = 0 (don't feed products)
    For inerts (not in reaction): ratio = 1.0
    
    Parameters:
        processes: Dictionary of process_id -> Process
        components: Dictionary of component_id -> Component
        
    Returns:
        Dictionary of component_id -> ratio relative to reference reactant
    """
    stoich_ratios = {}
    reference_comp = None
    reference_stoich = None
    
    # Find reactors and build stoichiometric ratios
    for proc in processes.values():
        if proc.process_type == ProcessType.REACTOR:
            for cid, coeff in proc.stoich.items():
                if coeff < 0:  # Reactant
                    stoich_ratios[cid] = abs(coeff)
                    # Use smallest stoich coefficient as reference
                    if reference_stoich is None or abs(coeff) < reference_stoich:
                        reference_comp = cid
                        reference_stoich = abs(coeff)
                elif coeff > 0:  # Product - don't feed
                    if cid not in stoich_ratios:
                        stoich_ratios[cid] = 0
    
    # Normalize ratios relative to reference
    if reference_comp and reference_stoich:
        for cid in stoich_ratios:
            if stoich_ratios[cid] > 0:
                stoich_ratios[cid] = stoich_ratios[cid] / reference_stoich
    
    # For components not in reaction (inerts), use ratio 1
    for cid in components:
        if cid not in stoich_ratios:
            stoich_ratios[cid] = 1.0
    
    return stoich_ratios


def reverse_solve(
    components: Dict[str, Component],
    processes: Dict[str, Process],
    tubes: Dict[str, Tube],
    target_output: Dict[str, float],
    alpha_eff: float
) -> Tuple[Optional[Dict[str, Dict[str, float]]], Optional[Dict[str, Dict[str, float]]], Optional[float], Dict[str, float], List[str]]:
    """
    Reverse solve: given desired output, calculate required input.
    
    Uses stoichiometric ratios to determine proper input proportions,
    then scales to achieve the target output.
    
    Parameters:
        components: Dictionary of component_id -> Component
        processes: Dictionary of process_id -> Process  
        tubes: Dictionary of tube_id -> Tube
        target_output: Dictionary of component_id -> target moles at output
        alpha_eff: Overall conversion efficiency
        
    Returns:
        results: Dictionary of tube_id -> {component_id -> moles} (or None if failed)
        required_input: Dictionary of process_id -> {component_id -> moles} (or None if failed)
        scale_factor: The scale factor used (or None if failed)
        stoich_ratios: Dictionary of component_id -> stoichiometric ratio
        warnings: List of warning messages
    """
    warnings = []
    
    # Find input and output processes
    input_procs = [p for p in processes.values() if p.process_type == ProcessType.INPUT]
    output_procs = [p for p in processes.values() if p.process_type == ProcessType.OUTPUT]
    
    if not input_procs:
        warnings.append("No INPUT processes found")
        return None, None, None, {}, warnings
    
    if not output_procs:
        warnings.append("No OUTPUT processes found")
        return None, None, None, {}, warnings
    
    # Get stoichiometric input ratios
    stoich_ratios = calculate_stoichiometric_input_ratios(processes, components)
    
    # Set stoichiometric input (1 mol of reference component equivalent)
    unit_input_moles = {}
    for inp in input_procs:
        unit_input_moles[inp.id] = {}
        for cid in components:
            unit_input_moles[inp.id][cid] = stoich_ratios.get(cid, 0)
    
    # Solve with stoichiometric input
    results, solve_warnings = solve_network(
        components,
        processes,
        tubes,
        unit_input_moles,
        alpha_eff
    )
    warnings.extend(solve_warnings)
    
    if not results:
        warnings.append("Could not solve network with unit input")
        return None, None, None, stoich_ratios, warnings
    
    # Check for negative values in unit solution
    for tid, comps in results.items():
        for cid, val in comps.items():
            if val < -1e-9:
                warnings.append("Base solution has negative flows - network may be infeasible")
                break
    
    # Find output tubes
    output_tube_ids = [tid for tid, tube in tubes.items() 
                       if tube.to_process in [p.id for p in output_procs]]
    
    # Calculate output per unit input
    unit_output = {}
    for cid in components:
        unit_output[cid] = sum(results.get(tid, {}).get(cid, 0) for tid in output_tube_ids)
    
    # Calculate scale factor - use the component with non-zero target and unit output
    scale_factors = {}
    for cid in components:
        if target_output.get(cid, 0) > 0 and abs(unit_output.get(cid, 0)) > 1e-9:
            scale_factors[cid] = target_output[cid] / unit_output[cid]
    
    if not scale_factors:
        warnings.append("Cannot determine scaling factor - no matching target/output components")
        return None, None, None, stoich_ratios, warnings
    
    # Use average scale factor
    avg_scale = sum(scale_factors.values()) / len(scale_factors)
    
    # Warn if scale factors are inconsistent
    if len(scale_factors) > 1:
        min_scale = min(scale_factors.values())
        max_scale = max(scale_factors.values())
        if max_scale > 0 and (max_scale - min_scale) / max_scale > 0.05:
            warnings.append(f"Target ratios don't match expected output ratios (scale factors vary: {min_scale:.4f} to {max_scale:.4f})")
    
    # Scale all results
    for tid in results:
        for cid in results[tid]:
            results[tid][cid] *= avg_scale
    
    # Calculate required input
    required_input = {}
    for pid in unit_input_moles:
        required_input[pid] = {}
        for cid in unit_input_moles[pid]:
            required_input[pid][cid] = unit_input_moles[pid][cid] * avg_scale
    
    # Check for negative results after scaling
    for tid, comps in results.items():
        for cid, val in comps.items():
            if val < -1e-9:
                warnings.append(f"Scaled solution has negative flows - target may be infeasible")
                break
    
    return results, required_input, avg_scale, stoich_ratios, warnings
