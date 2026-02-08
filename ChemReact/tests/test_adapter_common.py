import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "adapters" / "adapter_common.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("adapter_common", str(MODULE_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError("failed to load adapter_common")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestAdapterCommon(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = _load_module()
        cls.tmp_root = Path(tempfile.mkdtemp(prefix="chemreact_adapter_tests_", dir=str(ROOT)))

    def test_validate_request_missing_required(self):
        with self.assertRaises(ValueError):
            self.mod._validate_request({"target_smiles": "CCO"})

    def test_validate_request_type_error(self):
        with self.assertRaises(ValueError):
            self.mod._validate_request({"target_smiles": "CCO", "routes_file": "r.json", "top_k": "1"})

    def test_resolve_path_relative(self):
        request_file = ROOT / "samples" / "host_request.json"
        resolved = self.mod._resolve_path(request_file, "routes.json")
        self.assertTrue(resolved.endswith(str(Path("samples") / "routes.json")))

    def test_build_command_with_defaults(self):
        run_script = ROOT / "scripts" / "run_closed_loop.py"
        request_file = ROOT / "samples" / "host_request.json"
        payload = {"target_smiles": "CCO", "routes_file": "routes.json"}
        cmd = self.mod._build_command(run_script, request_file, payload, "outputs/default_run")
        self.assertIn("--target-smiles", cmd)
        self.assertIn("--routes-file", cmd)
        self.assertIn("--output-dir", cmd)
        out_idx = cmd.index("--output-dir") + 1
        self.assertIn("outputs", cmd[out_idx])

    def test_validate_request_auto_propose_without_routes(self):
        self.mod._validate_request({"target_smiles": "CCO", "auto_propose_routes": True})

    def test_validate_request_rejects_non_host_backend(self):
        with self.assertRaises(ValueError):
            self.mod._validate_request(
                {"target_smiles": "CCO", "auto_propose_routes": True, "planner_backend": "heuristic"}
            )

    def test_validate_request_rejects_auto_propose_max_routes_out_of_range(self):
        with self.assertRaises(ValueError):
            self.mod._validate_request(
                {"target_smiles": "CCO", "auto_propose_routes": True, "auto_propose_max_routes": 8}
            )

    def test_validate_request_rejects_planner_max_iterations_out_of_range(self):
        with self.assertRaises(ValueError):
            self.mod._validate_request(
                {"target_smiles": "CCO", "auto_propose_routes": True, "planner_max_iterations": 8}
            )

    def test_validate_request_rejects_fallback_host_command(self):
        with self.assertRaises(ValueError):
            self.mod._validate_request(
                {
                    "target_smiles": "CCO",
                    "auto_propose_routes": True,
                    "host_planner_command": "python skills/ChemReact/scripts/host_planner_fallback.py --prompt-file {prompt_file} --output-file {output_file} --request-file {request_file} --stage {stage} --round {round}",
                }
            )

    def test_validate_request_requires_full_host_placeholders(self):
        with self.assertRaises(ValueError):
            self.mod._validate_request(
                {
                    "target_smiles": "CCO",
                    "auto_propose_routes": True,
                    "host_planner_command": "my_host --prompt {prompt_file} --out {output_file}",
                }
            )

    def test_build_command_with_optional_fields(self):
        run_script = ROOT / "scripts" / "run_closed_loop.py"
        request_file = self.tmp_root / "request.json"
        request_file.write_text("{}", encoding="utf-8")
        payload = {
            "target_smiles": "CCO",
            "routes_file": str(ROOT / "samples" / "routes.json"),
            "strategy_file": str(ROOT / "samples" / "strategy.json"),
            "vis_plan_file": str(ROOT / "samples" / "vis_plan.json"),
            "top_k": 2,
            "target_legend": "Target",
            "force_field": "UFF",
            "emit_vis_plan_template": "tpl.json",
        }
        cmd = self.mod._build_command(run_script, request_file, payload, "outputs/default_run")
        self.assertIn("--strategy-file", cmd)
        self.assertIn("--vis-plan-file", cmd)
        self.assertIn("--top-k", cmd)
        self.assertIn("--force-field", cmd)
        self.assertIn("--emit-vis-plan-template", cmd)

    def test_build_command_auto_propose_host(self):
        run_script = ROOT / "scripts" / "run_closed_loop.py"
        request_file = ROOT / "samples" / "host_request.json"
        payload = {
            "target_smiles": "CCO",
            "auto_propose_routes": True,
            "planner_backend": "host",
            "auto_propose_max_routes": 5,
            "planner_request_only": True,
        }
        cmd = self.mod._build_command(run_script, request_file, payload, "outputs/default_run")
        self.assertIn("--auto-propose-routes", cmd)
        self.assertIn("--planner-backend", cmd)
        self.assertIn("--auto-propose-max-routes", cmd)
        self.assertIn("--planner-request-only", cmd)
        self.assertNotIn("--routes-file", cmd)

    def test_build_command_with_planner_repair_flags(self):
        run_script = ROOT / "scripts" / "run_closed_loop.py"
        request_file = ROOT / "samples" / "host_request.json"
        payload = {
            "target_smiles": "CCO",
            "auto_propose_routes": True,
            "planner_backend": "host",
            "planner_auto_repair": True,
            "planner_max_iterations": 4,
            "strict_audit": True,
            "host_planner_command": "my_host --prompt {prompt_file} --out {output_file}",
        }
        cmd = self.mod._build_command(run_script, request_file, payload, "outputs/default_run")
        self.assertIn("--planner-auto-repair", cmd)
        self.assertIn("--planner-max-iterations", cmd)
        self.assertIn("--strict-audit", cmd)
        self.assertIn("--host-planner-command", cmd)


if __name__ == "__main__":
    unittest.main()
