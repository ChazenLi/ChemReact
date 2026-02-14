"""
Chemical computation tools for MultiRetro.

Provides modular, composable tools for RDKit-based molecular operations:

- mol_parser:          SMILES parsing, normalization, and validation with LRU caching
- structure_analyzer:  Molecular structure analysis (features, atom-bond maps)
- functional_groups:   Functional group detection, protecting groups, danger alerts
- sa_scorer:           SA Score estimation and classification via SAThresholds
- atom_mapper:         RXNMapper wrapper for atom mapping reactions
- reaction_validator:  Reaction validation with 18 byproduct patterns, 3-tier verdict
- reaction_classifier: SMARTS + keyword reaction classification (33 rules from rules.json)
"""
