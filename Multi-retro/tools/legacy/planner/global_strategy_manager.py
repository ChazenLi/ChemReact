"""
MultiRetro Core Planner - Global Strategy Manager
==================================================
Maintains consistency across the entire retrosynthesis tree.
Prevents redundant protecting group cycles and enforces strategic constraints.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple
import json


# ─────────────────────────────────────────────────────────────────────────────
# Protecting Group Registry
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_PROTECTING_GROUPS: Dict[str, Dict[str, Any]] = {
    # Amine protecting groups
    "Boc": {
        "full_name": "tert-Butyloxycarbonyl",
        "target_fg": "amine",
        "install_reagents": ["Boc2O", "Et3N"],
        "remove_reagents": ["TFA", "HCl/dioxane"],
        "smarts_protected": "[N:1][C:2](=O)OC(C)(C)C",
        "orthogonal_to": ["Fmoc", "Cbz", "Alloc"],
    },
    "Fmoc": {
        "full_name": "9-Fluorenylmethyloxycarbonyl",
        "target_fg": "amine",
        "install_reagents": ["Fmoc-Cl", "NaHCO3"],
        "remove_reagents": ["piperidine/DMF"],
        "smarts_protected": "[N:1][C:2](=O)OCC1c2ccccc2c3ccccc13",
        "orthogonal_to": ["Boc", "Cbz", "Alloc"],
    },
    "Cbz": {
        "full_name": "Benzyloxycarbonyl",
        "target_fg": "amine",
        "install_reagents": ["Cbz-Cl", "NaOH"],
        "remove_reagents": ["H2/Pd-C", "HBr/AcOH"],
        "smarts_protected": "[N:1][C:2](=O)OCc1ccccc1",
        "orthogonal_to": ["Boc", "Fmoc", "TBS"],
    },
    # Alcohol protecting groups
    "TBS": {
        "full_name": "tert-Butyldimethylsilyl",
        "target_fg": "alcohol",
        "install_reagents": ["TBSCl", "imidazole"],
        "remove_reagents": ["TBAF", "HF/pyridine"],
        "smarts_protected": "[O:1][Si](C)(C)C(C)(C)C",
        "orthogonal_to": ["Boc", "Fmoc", "Cbz"],
    },
    "TIPS": {
        "full_name": "Triisopropylsilyl",
        "target_fg": "alcohol",
        "install_reagents": ["TIPSCl", "imidazole"],
        "remove_reagents": ["TBAF", "HF/pyridine"],
        "smarts_protected": "[O:1][Si](C(C)C)(C(C)C)C(C)C",
        "orthogonal_to": ["Boc", "Fmoc"],
    },
    "Bn": {
        "full_name": "Benzyl",
        "target_fg": "alcohol",
        "install_reagents": ["BnBr", "NaH"],
        "remove_reagents": ["H2/Pd-C"],
        "smarts_protected": "[O:1]Cc1ccccc1",
        "orthogonal_to": ["TBS", "TIPS"],
    },
    # Carboxylic acid protecting groups
    "tBu_ester": {
        "full_name": "tert-Butyl ester",
        "target_fg": "carboxylic_acid",
        "install_reagents": ["isobutene", "H2SO4"],
        "remove_reagents": ["TFA"],
        "smarts_protected": "[C:1](=O)OC(C)(C)C",
        "orthogonal_to": ["Bn_ester", "Me_ester"],
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Preferred Reaction Methods
# ─────────────────────────────────────────────────────────────────────────────

PREFERRED_METHODS: Dict[str, List[str]] = {
    "C-C_aryl-aryl": ["Suzuki", "Negishi", "Kumada"],
    "C-C_aryl-alkyl": ["Negishi", "Kumada", "Suzuki"],
    "C-C_alkyl-alkyl": ["Alkylation", "Grignard", "Aldol"],
    "C-N_aryl-N": ["Buchwald-Hartwig", "Ullmann", "Chan-Lam"],
    "C-N_alkyl-N": ["Reductive_amination", "N-alkylation", "Gabriel"],
    "C-O_ester": ["Fischer_esterification", "Steglich", "Yamaguchi"],
    "C-O_ether": ["Williamson", "Mitsunobu"],
    "amide": ["EDC_coupling", "DCC_coupling", "Schotten-Baumann"],
}


# ─────────────────────────────────────────────────────────────────────────────
# Protecting Group Operation
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class PGOperation:
    """Record of a protecting group operation."""
    step_index: int
    node_id: str
    pg_type: str
    operation: str  # "install" or "remove"
    target_fg: str
    smiles_before: str
    smiles_after: str
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "step_index": self.step_index,
            "node_id": self.node_id,
            "pg_type": self.pg_type,
            "operation": self.operation,
            "target_fg": self.target_fg,
            "smiles_before": self.smiles_before,
            "smiles_after": self.smiles_after,
        }


# ─────────────────────────────────────────────────────────────────────────────
# Global Strategy Manager
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class GlobalStrategyConfig:
    """Configuration for global strategy constraints."""
    max_pg_operations: int = 2
    no_redundant_cycles: bool = True
    allowed_pg_types: List[str] = field(default_factory=lambda: ["Boc", "Fmoc", "TBS", "TIPS"])
    preferred_methods: Dict[str, List[str]] = field(default_factory=dict)
    forbidden_reagents: List[str] = field(default_factory=list)
    max_steps: int = 10
    convergent_preferred: bool = True


class GlobalStrategyManager:
    """
    Maintains consistency across the entire retrosynthesis tree.
    
    Key responsibilities:
    1. Track protecting group operations
    2. Detect redundant protection/deprotection cycles
    3. Enforce preferred reaction methods
    4. Propagate constraints to child nodes
    """
    
    def __init__(self, config: Optional[GlobalStrategyConfig] = None):
        self.config = config or GlobalStrategyConfig()
        self.pg_history: List[PGOperation] = []
        self.step_count = 0
        self.issues: List[Dict[str, Any]] = []
        self._active_pg: Dict[str, Set[str]] = {}  # node_id -> set of active PGs
    
    def get_strategy_context(self) -> Dict[str, Any]:
        """Get current strategy context for prompts."""
        return {
            "protecting_group_policy": {
                "allowed_groups": self.config.allowed_pg_types,
                "max_pg_operations": self.config.max_pg_operations,
                "current_pg_count": len(self.pg_history),
                "no_redundant_cycles": self.config.no_redundant_cycles,
            },
            "preferred_methods": self.config.preferred_methods or PREFERRED_METHODS,
            "forbidden_reagents": self.config.forbidden_reagents,
            "current_step_count": self.step_count,
            "max_steps": self.config.max_steps,
            "convergent_preferred": self.config.convergent_preferred,
            "active_protecting_groups": {k: list(v) for k, v in self._active_pg.items()},
        }
    
    def register_pg_operation(self, operation: PGOperation) -> List[str]:
        """
        Register a protecting group operation and check for issues.
        
        Returns:
            List of issue codes if problems detected
        """
        issues = []
        
        # Check if PG type is allowed
        if operation.pg_type not in self.config.allowed_pg_types:
            issues.append(f"FORBIDDEN_PG_TYPE:{operation.pg_type}")
        
        # Check max PG operations
        if len(self.pg_history) >= self.config.max_pg_operations:
            issues.append(f"MAX_PG_OPERATIONS_EXCEEDED:{len(self.pg_history)+1}/{self.config.max_pg_operations}")
        
        # Check for redundant cycles
        if self.config.no_redundant_cycles:
            cycle_issue = self._check_redundant_cycle(operation)
            if cycle_issue:
                issues.append(cycle_issue)
        
        # Register the operation
        self.pg_history.append(operation)
        
        # Update active PG tracking
        if operation.operation == "install":
            if operation.node_id not in self._active_pg:
                self._active_pg[operation.node_id] = set()
            self._active_pg[operation.node_id].add(operation.pg_type)
        elif operation.operation == "remove":
            if operation.node_id in self._active_pg:
                self._active_pg[operation.node_id].discard(operation.pg_type)
        
        # Log issues
        for issue in issues:
            self.issues.append({
                "code": issue,
                "step_index": operation.step_index,
                "pg_type": operation.pg_type,
                "operation": operation.operation,
            })
        
        return issues
    
    def _check_redundant_cycle(self, new_op: PGOperation) -> Optional[str]:
        """
        Check if this operation creates a redundant cycle.
        
        Redundant cycle example:
        - Step 2: Install Boc on NodeA
        - Step 4: Remove Boc from NodeA
        - Step 6: Install Boc on NodeA (REDUNDANT!)
        """
        if new_op.operation != "install":
            return None
        
        # Look for previous remove of same PG type on same functional group
        for prev_op in reversed(self.pg_history):
            if (prev_op.operation == "remove" and 
                prev_op.pg_type == new_op.pg_type and
                prev_op.target_fg == new_op.target_fg):
                
                # Found a deprotection, now look for earlier protection
                for earlier_op in self.pg_history:
                    if earlier_op.step_index >= prev_op.step_index:
                        break
                    if (earlier_op.operation == "install" and
                        earlier_op.pg_type == new_op.pg_type and
                        earlier_op.target_fg == new_op.target_fg):
                        return f"REDUNDANT_PG_CYCLE:{new_op.pg_type}"
        
        return None
    
    def validate_step(self, step: Dict[str, Any]) -> List[str]:
        """
        Validate a proposed step against global strategy.
        
        Args:
            step: Step dict with reaction_class, reagents, etc.
            
        Returns:
            List of issue codes
        """
        issues = []
        
        # Check forbidden reagents
        reagents = step.get("reagents", [])
        for reagent in reagents:
            if isinstance(reagent, str):
                for forbidden in self.config.forbidden_reagents:
                    if forbidden.lower() in reagent.lower():
                        issues.append(f"FORBIDDEN_REAGENT:{reagent}")
        
        # Check max steps
        if self.step_count >= self.config.max_steps:
            issues.append(f"MAX_STEPS_EXCEEDED:{self.step_count+1}/{self.config.max_steps}")
        
        # Check preferred methods
        reaction_class = step.get("reaction_class", "").lower()
        bond_type = self._infer_bond_type(step)
        if bond_type and bond_type in PREFERRED_METHODS:
            preferred = PREFERRED_METHODS[bond_type]
            is_preferred = any(p.lower() in reaction_class for p in preferred)
            if not is_preferred:
                issues.append(f"NON_PREFERRED_METHOD:{reaction_class}")
        
        self.step_count += 1
        return issues
    
    def _infer_bond_type(self, step: Dict[str, Any]) -> Optional[str]:
        """Infer bond type from step information."""
        reaction_class = step.get("reaction_class", "").lower()
        
        if any(x in reaction_class for x in ["suzuki", "negishi", "kumada", "heck"]):
            return "C-C_aryl-aryl"
        if any(x in reaction_class for x in ["buchwald", "ullmann", "amination"]):
            return "C-N_aryl-N"
        if any(x in reaction_class for x in ["amide", "coupling", "edc", "dcc"]):
            return "amide"
        if any(x in reaction_class for x in ["ester", "fischer"]):
            return "C-O_ester"
        if any(x in reaction_class for x in ["williamson", "ether"]):
            return "C-O_ether"
        
        return None
    
    def get_pg_recommendations(self, functional_groups: List[str]) -> List[Dict[str, Any]]:
        """
        Get protecting group recommendations for functional groups.
        
        Args:
            functional_groups: List of FG types needing protection
            
        Returns:
            List of PG recommendations with orthogonality info
        """
        recommendations = []
        
        for fg in functional_groups:
            candidates = []
            for pg_name, pg_info in KNOWN_PROTECTING_GROUPS.items():
                if pg_info["target_fg"] == fg and pg_name in self.config.allowed_pg_types:
                    # Check if orthogonal to currently active PGs
                    active_pgs = set()
                    for node_pgs in self._active_pg.values():
                        active_pgs.update(node_pgs)
                    
                    is_orthogonal = all(
                        pg_name in KNOWN_PROTECTING_GROUPS.get(active, {}).get("orthogonal_to", [])
                        for active in active_pgs
                    )
                    
                    candidates.append({
                        "pg_type": pg_name,
                        "full_name": pg_info["full_name"],
                        "install_reagents": pg_info["install_reagents"],
                        "remove_reagents": pg_info["remove_reagents"],
                        "is_orthogonal": is_orthogonal,
                        "priority": 1 if is_orthogonal else 2,
                    })
            
            # Sort by orthogonality
            candidates.sort(key=lambda x: x["priority"])
            
            if candidates:
                recommendations.append({
                    "functional_group": fg,
                    "recommended_pgs": candidates[:2],  # Top 2
                })
        
        return recommendations
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary of strategy state."""
        return {
            "total_pg_operations": len(self.pg_history),
            "pg_history": [op.to_dict() for op in self.pg_history],
            "total_steps": self.step_count,
            "issues": self.issues,
            "active_protecting_groups": {k: list(v) for k, v in self._active_pg.items()},
            "config": {
                "max_pg_operations": self.config.max_pg_operations,
                "allowed_pg_types": self.config.allowed_pg_types,
                "no_redundant_cycles": self.config.no_redundant_cycles,
                "max_steps": self.config.max_steps,
            },
        }


# ─────────────────────────────────────────────────────────────────────────────
# Factory Functions
# ─────────────────────────────────────────────────────────────────────────────

def create_strategy_manager(
    constraints: Optional[Dict[str, Any]] = None
) -> GlobalStrategyManager:
    """
    Create a GlobalStrategyManager with optional constraints.
    
    Args:
        constraints: Optional dict with config overrides
        
    Returns:
        Configured GlobalStrategyManager
    """
    config = GlobalStrategyConfig()
    
    if constraints:
        if "max_pg_operations" in constraints:
            config.max_pg_operations = constraints["max_pg_operations"]
        if "allowed_pg_types" in constraints:
            config.allowed_pg_types = constraints["allowed_pg_types"]
        if "no_redundant_cycles" in constraints:
            config.no_redundant_cycles = constraints["no_redundant_cycles"]
        if "preferred_methods" in constraints:
            config.preferred_methods = constraints["preferred_methods"]
        if "forbidden_reagents" in constraints:
            config.forbidden_reagents = constraints["forbidden_reagents"]
        if "max_steps" in constraints:
            config.max_steps = constraints["max_steps"]
        if "convergent_preferred" in constraints:
            config.convergent_preferred = constraints["convergent_preferred"]
    
    return GlobalStrategyManager(config)


# ─────────────────────────────────────────────────────────────────────────────
# CLI / Testing Entry Point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Demo usage
    manager = create_strategy_manager({
        "max_pg_operations": 2,
        "allowed_pg_types": ["Boc", "Fmoc", "TBS"],
    })
    
    # Simulate operations
    op1 = PGOperation(
        step_index=1,
        node_id="node_001",
        pg_type="Boc",
        operation="install",
        target_fg="amine",
        smiles_before="Nc1ccccc1",
        smiles_after="CC(C)(C)OC(=O)Nc1ccccc1",
    )
    
    issues1 = manager.register_pg_operation(op1)
    print(f"Op1 issues: {issues1}")
    
    # Get recommendations
    recs = manager.get_pg_recommendations(["alcohol"])
    print(f"PG recommendations: {json.dumps(recs, indent=2)}")
    
    # Get summary
    print(f"Summary: {json.dumps(manager.get_summary(), indent=2)}")
