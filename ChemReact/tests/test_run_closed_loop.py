import json
import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_closed_loop.py"
SAMPLES = ROOT / "samples"


def _load_run_closed_loop_module():
    spec = importlib.util.spec_from_file_location("run_closed_loop_mod", str(SCRIPT))
    if spec is None or spec.loader is None:
        raise RuntimeError("failed to load run_closed_loop.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestRunClosedLoopModes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp_root = Path(tempfile.mkdtemp(prefix="chemreact_tests_", dir=str(ROOT)))
        cls.mod = _load_run_closed_loop_module()

    def _run(self, args):
        proc = subprocess.run(
            [sys.executable, str(SCRIPT)] + args,
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, msg=proc.stdout + "\n" + proc.stderr)
        return proc

    def _run_raw(self, args):
        return subprocess.run(
            [sys.executable, str(SCRIPT)] + args,
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=False,
        )

    def test_template_only_mode(self):
        out = self.tmp_root / "template"
        self._run(["--output-dir", str(out), "--emit-vis-plan-template-only"])
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertEqual(summary["mode"], "template_only")
        self.assertIn("timings", summary)
        self.assertTrue((out / "schemas" / "routes.schema.json").exists())

    def test_validate_only_mode(self):
        out = self.tmp_root / "validate"
        self._run(
            [
                "--target-smiles",
                "CCO",
                "--routes-file",
                str(SAMPLES / "routes.json"),
                "--strategy-file",
                str(SAMPLES / "strategy.json"),
                "--vis-plan-file",
                str(SAMPLES / "vis_plan.json"),
                "--output-dir",
                str(out),
                "--validate-only",
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertEqual(summary["mode"], "validate_only")
        self.assertIn("timings", summary)

    def test_full_mode(self):
        out = self.tmp_root / "full"
        self._run(
            [
                "--target-smiles",
                "CCO",
                "--routes-file",
                str(SAMPLES / "routes.json"),
                "--strategy-file",
                str(SAMPLES / "strategy.json"),
                "--vis-plan-file",
                str(SAMPLES / "vis_plan.json"),
                "--output-dir",
                str(out),
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertIn("timings", summary)
        self.assertIn("stages", summary["timings"])
        self.assertTrue((out / "RETRO_REPORT.md").exists())
        self.assertTrue((out / "images" / "target.png").exists())
        self.assertIn("functional_groups", summary["rdkit"])

    def test_planner_request_only_mode(self):
        out = self.tmp_root / "planner_req"
        self._run(
            [
                "--target-smiles",
                "CCO",
                "--planner-request-only",
                "--planner-backend",
                "host",
                "--output-dir",
                str(out),
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertEqual(summary["mode"], "planner_request_only")
        self.assertTrue((out / "planner_request.json").exists())
        self.assertTrue((out / "planner_routes.template.json").exists())
        self.assertTrue((out / "planner_prompt.txt").exists())
        planner_request = json.loads((out / "planner_request.json").read_text(encoding="utf-8"))
        self.assertIn("target_chemistry", planner_request)
        self.assertIn("functional_groups", planner_request["target_chemistry"])
        self.assertIn("ring_profile", planner_request["target_chemistry"])

    def test_auto_propose_min_routes_is_three(self):
        out = self.tmp_root / "planner_req_auto_min3"
        self._run(
            [
                "--target-smiles",
                "CCO",
                "--planner-request-only",
                "--auto-propose-routes",
                "--auto-propose-max-routes",
                "3",
                "--output-dir",
                str(out),
            ]
        )
        planner_request = json.loads((out / "planner_request.json").read_text(encoding="utf-8"))
        self.assertEqual(planner_request["max_routes"], 3)

    def test_auto_propose_max_routes_is_clamped_to_five(self):
        out = self.tmp_root / "planner_req_auto_max5"
        self._run(
            [
                "--target-smiles",
                "CCO",
                "--planner-request-only",
                "--auto-propose-routes",
                "--auto-propose-max-routes",
                "9",
                "--output-dir",
                str(out),
            ]
        )
        planner_request = json.loads((out / "planner_request.json").read_text(encoding="utf-8"))
        self.assertEqual(planner_request["max_routes"], 5)

    def test_auto_propose_host_single_command_closed_loop(self):
        out = self.tmp_root / "auto_host_closed_loop"
        out.mkdir(parents=True, exist_ok=True)
        planner_routes = out / "planner_routes.json"
        shutil.copyfile(SAMPLES / "routes.json", planner_routes)
        self._run(
            [
                "--target-smiles",
                "CCO",
                "--auto-propose-routes",
                "--planner-backend",
                "host",
                "--no-planner-auto-repair",
                "--output-dir",
                str(out),
                "--validate-only",
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertEqual(summary["mode"], "validate_only")
        self.assertEqual(Path(summary["validated"]["routes_file"]).resolve(), planner_routes.resolve())

    def test_auto_propose_defaults_to_auto_repair(self):
        out = self.tmp_root / "auto_repair_default"
        out.mkdir(parents=True, exist_ok=True)
        planner_routes = out / "planner_routes.json"
        planner_routes_payload = [
            {
                "route_id": 1,
                "score": 8.0,
                "audit_verdict": "PASS",
                "critical_issues": [],
                "precursors": ["CCO"],
                "steps": [
                    {
                        "reaction_class": "Bad Step",
                        "reagents": ["X"],
                        "conditions": "RT",
                        "reaction_smiles": "CCO>>CCN"
                    }
                ]
            }
        ]
        planner_routes.write_text(json.dumps(planner_routes_payload, ensure_ascii=False, indent=2), encoding="utf-8")
        self._run(
            [
                "--target-smiles",
                "CCO",
                "--auto-propose-routes",
                "--planner-backend",
                "host",
                "--planner-max-iterations",
                "2",
                "--output-dir",
                str(out),
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertEqual(summary["mode"], "planner_repair_requested")
        self.assertIn("planner_repair_prompt", summary["outputs"])

    def test_auto_propose_rejects_heuristic_backend(self):
        out = self.tmp_root / "invalid_backend"
        proc = self._run_raw(
            [
                "--target-smiles",
                "CCO",
                "--auto-propose-routes",
                "--planner-backend",
                "heuristic",
                "--output-dir",
                str(out),
                "--validate-only",
            ]
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("invalid choice", proc.stdout + proc.stderr)

    def test_strict_audit_blocks_invalid_pass_route(self):
        out = self.tmp_root / "strict_gate"
        out.mkdir(parents=True, exist_ok=True)
        routes_path = out / "routes_invalid.json"
        routes_payload = [
            {
                "route_id": 1,
                "score": 9.2,
                "audit_verdict": "PASS",
                "critical_issues": [],
                "precursors": ["CCO"],
                "steps": [
                    {
                        "reaction_class": "Invalid Mapping",
                        "reagents": ["X"],
                        "conditions": "RT",
                        "reaction_smiles": "CCO>>CCN"
                    }
                ]
            }
        ]
        routes_path.write_text(json.dumps(routes_payload, ensure_ascii=False, indent=2), encoding="utf-8")

        self._run(
            [
                "--target-smiles",
                "CCO",
                "--routes-file",
                str(routes_path),
                "--strategy-file",
                str(SAMPLES / "strategy.json"),
                "--output-dir",
                str(out),
                "--strict-audit",
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertEqual(summary["selected_route_ids"], [])
        self.assertEqual(summary["rejected_route_ids"], [1])
        self.assertEqual(summary["audit"]["strict_mode"], True)
        self.assertEqual(summary["audit"]["rejected_routes"], 1)
        self.assertEqual(summary["audit"]["analysis_required_for_pass"], True)
        report_text = (out / "RETRO_REPORT.md").read_text(encoding="utf-8")
        self.assertIn("Quality Score", report_text)
        self.assertTrue((out / "images" / "route_1_step_1.png").exists())

    def test_strict_audit_hydrogen_conservation_is_warned_for_small_molecule_loss(self):
        out = self.tmp_root / "strict_hydrogen"
        out.mkdir(parents=True, exist_ok=True)
        routes_path = out / "routes_hydrogen_invalid.json"
        routes_payload = [
            {
                "route_id": 1,
                "score": 8.0,
                "audit_verdict": "PASS",
                "critical_issues": [],
                "precursors": ["CCO"],
                "steps": [
                    {
                        "reaction_class": "Oxidation",
                        "reagents": ["oxidant"],
                        "conditions": "RT",
                        "reaction_smiles": "CCO>>CC=O"
                    }
                ]
            }
        ]
        routes_path.write_text(json.dumps(routes_payload, ensure_ascii=False, indent=2), encoding="utf-8")

        self._run(
            [
                "--target-smiles",
                "CC=O",
                "--routes-file",
                str(routes_path),
                "--output-dir",
                str(out),
                "--strict-audit",
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["selected_route_ids"], [1])
        self.assertEqual(summary["rejected_route_ids"], [])
        report_text = (out / "RETRO_REPORT.md").read_text(encoding="utf-8")
        self.assertIn("ATOM_BALANCE_INFERRED_BYPRODUCT", report_text)

    def test_strict_audit_can_infer_water_byproduct_for_esterification(self):
        out = self.tmp_root / "strict_infer_byproduct"
        out.mkdir(parents=True, exist_ok=True)
        routes_path = out / "routes_infer_byproduct.json"
        routes_payload = [
            {
                "route_id": 1,
                "score": 8.0,
                "audit_verdict": "PASS",
                "critical_issues": [],
                "precursors": ["CC(=O)O", "CO"],
                "steps": [
                    {
                        "reaction_class": "Esterification",
                        "reagents": ["acid catalyst"],
                        "conditions": "reflux",
                        "reaction_smiles": "CC(=O)O.CO>>CC(=O)OC"
                    }
                ]
            }
        ]
        routes_path.write_text(json.dumps(routes_payload, ensure_ascii=False, indent=2), encoding="utf-8")

        self._run(
            [
                "--target-smiles",
                "CC(=O)OC",
                "--routes-file",
                str(routes_path),
                "--output-dir",
                str(out),
                "--strict-audit",
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["selected_route_ids"], [1])
        report_text = (out / "RETRO_REPORT.md").read_text(encoding="utf-8")
        self.assertIn("inferred byproducts", report_text.lower())

    def test_planner_auto_repair_requests_next_round_when_no_host_command(self):
        out = self.tmp_root / "planner_repair_requested"
        out.mkdir(parents=True, exist_ok=True)
        planner_routes = out / "planner_routes.json"
        planner_routes_payload = [
            {
                "route_id": 1,
                "score": 8.0,
                "audit_verdict": "PASS",
                "critical_issues": [],
                "precursors": ["CCO"],
                "steps": [
                    {
                        "reaction_class": "Bad Step",
                        "reagents": ["X"],
                        "conditions": "RT",
                        "reaction_smiles": "CCO>>CCN"
                    }
                ]
            }
        ]
        planner_routes.write_text(json.dumps(planner_routes_payload, ensure_ascii=False, indent=2), encoding="utf-8")
        self._run(
            [
                "--target-smiles",
                "CCO",
                "--auto-propose-routes",
                "--planner-backend",
                "host",
                "--planner-auto-repair",
                "--planner-max-iterations",
                "2",
                "--output-dir",
                str(out),
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertEqual(summary["mode"], "planner_repair_requested")
        self.assertIn("planner_repair_prompt", summary["outputs"])
        self.assertIn("planner_loop", summary["outputs"])

    def test_auto_propose_strict_enforces_style_mix(self):
        out = self.tmp_root / "style_mix_gate"
        out.mkdir(parents=True, exist_ok=True)
        planner_routes = out / "planner_routes.json"
        planner_routes_payload = []
        for rid in [1, 2, 3, 4]:
            planner_routes_payload.append(
                {
                    "route_id": rid,
                    "synthesis_style": "convergent",
                    "score": 7.0,
                    "audit_verdict": "PASS",
                    "critical_issues": [],
                    "precursors": ["CCBr", "[OH-]"],
                    "steps": [
                        {
                            "reaction_class": "Substitution",
                            "reagents": ["NaOH"],
                            "conditions": "RT",
                            "reaction_smiles": "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][OH:4].[Br-:3]"
                        }
                    ],
                }
            )
        planner_routes.write_text(json.dumps(planner_routes_payload, ensure_ascii=False, indent=2), encoding="utf-8")

        self._run(
            [
                "--target-smiles",
                "CCO",
                "--auto-propose-routes",
                "--planner-backend",
                "host",
                "--no-planner-auto-repair",
                "--output-dir",
                str(out),
                "--validate-only",
                "--strict-audit",
            ]
        )
        summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["success"])
        self.assertEqual(summary["mode"], "validate_only")
        self.assertEqual(summary["audit"]["portfolio_style_issue"] is not None, True)
        self.assertEqual(
            summary["audit"]["portfolio_style_requirement"],
            "at least 3 routes with mixed convergent/linear styles",
        )

    def test_step_audit_has_structured_reaction_type_and_atom_balance(self):
        class AtomUtilsStub:
            @staticmethod
            def validate_mapped_reaction(_):
                return {"valid": True, "completeness": 1.0}

            @staticmethod
            def get_bond_changes(_):
                return {"formed": [], "broken": []}

        step = {
            "reaction_class": "Esterification",
            "reaction_smiles": "CC(=O)O.CO>>CC(=O)OC",
        }
        result = self.mod._audit_single_step(
            step=step,
            step_index=1,
            route_precursors_norm=[],
            mapper=None,
            atom_utils_mod=AtomUtilsStub(),
        )
        self.assertIn("reaction_type", result)
        self.assertIn("atom_balance", result)
        self.assertEqual(result["reaction_type"]["id"], "esterification")
        self.assertIn("byproduct_candidates", result["reaction_type"])
        self.assertGreater(result["reaction_type"]["byproduct_candidates"][0]["confidence"], 0.0)
        self.assertIn("O", result["atom_balance"]["inferred_byproducts"])
        self.assertIn("inferred_byproduct_details", result["atom_balance"])
        self.assertGreater(result["atom_balance"]["inferred_byproduct_details"][0]["confidence"], 0.0)
        self.assertGreater(result["atom_balance"]["inferred_byproduct_confidence"], 0.0)
        self.assertEqual(result["atom_balance"]["balanced"], True)

    def test_repair_prompt_includes_step_level_reaction_type_guidance(self):
        repair_request = {
            "target_smiles": "CCO",
            "rejected_routes": [
                {
                    "route_id": 7,
                    "audit_details": {
                        "hard_issues": ["Atom conservation failed between reactants and products"],
                        "soft_issues": [],
                        "steps": [
                            {
                                "step_index": 1,
                                "reaction_type": {
                                    "id": "esterification",
                                    "confidence": 0.82,
                                    "source": "smarts+keyword",
                                    "byproduct_candidates": [{"smiles": "O", "confidence": 0.77}],
                                },
                                "mapping": {
                                    "completeness": 0.61,
                                    "notes": ["Atom mapping completeness is low (0.61)"],
                                },
                                "atom_balance": {
                                    "missing_product_delta": {"H": 2, "O": 1},
                                    "excess_product_delta": {},
                                },
                            }
                        ],
                    },
                }
            ],
        }
        prompt = self.mod._build_repair_prompt(repair_request)
        self.assertIn("Step-level Reaction-Type Guidance", prompt)
        self.assertIn("Step 1: reaction_type=esterification", prompt)
        self.assertIn("suggested explicit byproducts", prompt)
        self.assertIn("output fully mapped rs>>ps", prompt)


if __name__ == "__main__":
    unittest.main()
